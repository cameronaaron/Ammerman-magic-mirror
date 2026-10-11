"""FastAPI application: upload a genome to an account, read your traits back.

This is the server half of the project's stated next step — "parse the 23andMe
files when they are uploaded to a server through a web front end … associated
with a specific user so that when the person logs into their account on the
mirror, their DNA info can be displayed to them."

Endpoints (all JSON unless noted):

    POST /api/register              {username, password}          -> {user_id}
    POST /api/login                 {username, password}          -> {token}
    POST /api/genomes               multipart file (Bearer auth)  -> {genome_id…}
    GET  /api/genomes               (Bearer auth)                 -> [genome…]
    GET  /api/genomes/{id}/traits   (Bearer auth)                 -> {traits…}
    GET  /                          minimal HTML upload page

The mirror pulls a logged-in user's traits from ``/api/genomes/{id}/traits``.
SNPedia is loaded once and shared across requests.
"""

from __future__ import annotations

import io
from pathlib import Path

from fastapi import Depends, FastAPI, File, Header, HTTPException, Query, UploadFile
from fastapi.responses import HTMLResponse, Response
from PIL import Image
from pydantic import BaseModel

from .. import snpedia
from ..app import (
    SHARE_SIZES,
    Reflection,
    build_resolver,
    compose_share_image,
)
from ..config import PROJECT_ROOT, Settings
from ..genome import load_genome
from ..traits import build_traits
from .storage import Storage, StorageError

_MAX_UPLOAD_BYTES = 64 * 1024 * 1024  # 23andMe exports are ~15 MB; cap the abuse
_MAX_PHOTO_BYTES = 25 * 1024 * 1024


class Credentials(BaseModel):
    username: str
    password: str


class TraitService:
    """Lazily-loaded SNPedia table + resolver, reused across requests."""

    def __init__(self, settings: Settings) -> None:
        self.settings = settings
        self._table = None
        self._resolver = build_resolver(settings)

    def _ensure_table(self):
        if self._table is None:
            self._table = snpedia.resolve_table(self.settings, log=lambda *_: None)
        return self._table

    def reflection_for(self, path: Path) -> Reflection:
        """Load a genome and rank its traits against the cached SNPedia table."""
        table = self._ensure_table()
        genome = load_genome(path, resolver=self._resolver)
        matched_all = build_traits(genome, table, min_magnitude=0.0, require_summary=True)
        notable = build_traits(
            genome,
            table,
            min_magnitude=self.settings.min_magnitude,
            require_summary=True,
            limit=self.settings.max_traits,
        )
        return Reflection(traits=notable, genome=genome, matched_total=len(matched_all))

    def traits_for(self, path: Path) -> list[dict]:
        return [
            {
                "rsid": t.rsid,
                "genotype": t.display_genotype,
                "magnitude": t.magnitude,
                "repute": t.repute,
                "category": t.category,
                "summary": t.summary,
                "url": t.snpedia_url,
            }
            for t in self.reflection_for(path).traits
        ]


def create_app(
    settings: Settings | None = None, data_dir: str | Path | None = None
) -> FastAPI:
    """Build the FastAPI app.  ``data_dir`` overrides where accounts are stored."""
    settings = settings or Settings()
    storage = Storage(data_dir or (PROJECT_ROOT / "web_data"))
    service = TraitService(settings)

    app = FastAPI(title="Genetic Reflexions", version="2.0.0")
    app.state.storage = storage
    app.state.service = service

    def current_user(authorization: str = Header(default="")) -> int:
        token = authorization.removeprefix("Bearer ").strip()
        user_id = storage.user_for_token(token)
        if user_id is None:
            raise HTTPException(status_code=401, detail="invalid or missing token")
        return user_id

    @app.post("/api/register", status_code=201)
    def register(creds: Credentials) -> dict:
        try:
            user_id = storage.create_user(creds.username, creds.password)
        except StorageError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        return {"user_id": user_id, "username": creds.username.strip().lower()}

    @app.post("/api/login")
    def login(creds: Credentials) -> dict:
        user_id = storage.authenticate(creds.username, creds.password)
        if user_id is None:
            raise HTTPException(status_code=401, detail="bad username or password")
        return {"token": storage.create_session(user_id), "user_id": user_id}

    @app.post("/api/genomes", status_code=201)
    async def upload_genome(
        file: UploadFile = File(...), user_id: int = Depends(current_user)
    ) -> dict:
        data = await file.read(_MAX_UPLOAD_BYTES + 1)
        if len(data) > _MAX_UPLOAD_BYTES:
            raise HTTPException(status_code=413, detail="genome file too large")
        try:
            record = storage.save_genome(user_id, file.filename or "genome.txt", data)
        except StorageError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        return {
            "genome_id": record.id,
            "filename": record.filename,
            "markers": record.markers,
        }

    @app.get("/api/genomes")
    def list_genomes(user_id: int = Depends(current_user)) -> dict:
        return {
            "genomes": [
                {"id": g.id, "filename": g.filename, "markers": g.markers}
                for g in storage.list_genomes(user_id)
            ]
        }

    @app.get("/api/genomes/{genome_id}/traits")
    def genome_traits(genome_id: int, user_id: int = Depends(current_user)) -> dict:
        record = storage.get_genome(user_id, genome_id)
        if record is None or not record.stored_path.exists():
            raise HTTPException(status_code=404, detail="genome not found")
        return {"filename": record.filename, "traits": service.traits_for(record.stored_path)}

    @app.post("/api/genomes/{genome_id}/share")
    async def share_image(
        genome_id: int,
        file: UploadFile = File(...),
        fmt: str = Query("portrait", description="square · portrait · story · wide"),
        user_id: int = Depends(current_user),
    ) -> Response:
        if fmt not in SHARE_SIZES:
            raise HTTPException(status_code=400, detail=f"unknown format '{fmt}'")
        record = storage.get_genome(user_id, genome_id)
        if record is None or not record.stored_path.exists():
            raise HTTPException(status_code=404, detail="genome not found")
        data = await file.read(_MAX_PHOTO_BYTES + 1)
        if len(data) > _MAX_PHOTO_BYTES:
            raise HTTPException(status_code=413, detail="photo too large")
        try:
            photo = Image.open(io.BytesIO(data))
            photo.load()
        except Exception as exc:  # noqa: BLE001 - any decode failure is a 400
            raise HTTPException(status_code=400, detail="unreadable image") from exc

        reflection = service.reflection_for(record.stored_path)
        image = compose_share_image(
            settings, reflection, photo, fmt=fmt, detect=True, log=lambda *_: None
        )
        buffer = io.BytesIO()
        image.save(buffer, format="PNG")
        return Response(content=buffer.getvalue(), media_type="image/png")

    @app.get("/", response_class=HTMLResponse)
    def index() -> str:
        return _INDEX_HTML

    return app


_INDEX_HTML = """<!doctype html>
<title>Genetic Reflexions — upload</title>
<style>
 body{font:16px system-ui;max-width:40rem;margin:3rem auto;padding:0 1rem;
      background:#14161c;color:#e7e9ee}
 h1{font-weight:600} input,button{font:inherit;padding:.5rem;margin:.25rem 0}
 fieldset{border:1px solid #333;border-radius:8px;margin:1rem 0}
 pre{background:#0e0f13;padding:1rem;border-radius:8px;overflow:auto}
 a{color:#8ab4ff}
</style>
<h1>Genetic Reflexions</h1>
<p>Scaffold front end. Register, log in, then upload a 23andMe raw file. The
mirror reads your traits back from <code>/api/genomes/&lt;id&gt;/traits</code>.</p>
<fieldset><legend>1 · Account</legend>
 <input id=u placeholder=username> <input id=p type=password placeholder=password>
 <button onclick=call('register')>Register</button>
 <button onclick=call('login')>Log in</button>
</fieldset>
<fieldset><legend>2 · Upload genome</legend>
 <input type=file id=f> <button onclick=upload()>Upload</button>
 <button onclick=list()>My genomes</button>
</fieldset>
<fieldset><legend>3 · Make a shareable image</legend>
 genome id <input id=gid size=4 placeholder=1>
 <select id=fmt>
  <option>portrait</option><option>square</option><option>story</option><option>wide</option>
 </select>
 photo <input type=file id=photo accept=image/*>
 <button onclick=share()>Create</button>
 <div><img id=preview alt="" style="max-width:100%;margin-top:1rem;border-radius:8px"></div>
</fieldset>
<pre id=out>—</pre>
<script>
let token='';
const out=document.getElementById('out');
const show=o=>out.textContent=JSON.stringify(o,null,2);
async function call(kind){
 const r=await fetch('/api/'+kind,{method:'POST',headers:{'Content-Type':'application/json'},
   body:JSON.stringify({username:u.value,password:p.value})});
 const j=await r.json(); if(j.token)token=j.token; show(j);
}
async function upload(){
 const fd=new FormData(); fd.append('file',f.files[0]);
 const r=await fetch('/api/genomes',{method:'POST',headers:{Authorization:'Bearer '+token},body:fd});
 show(await r.json());
}
async function list(){
 const r=await fetch('/api/genomes',{headers:{Authorization:'Bearer '+token}});
 show(await r.json());
}
async function share(){
 const fd=new FormData(); fd.append('file',photo.files[0]);
 const r=await fetch('/api/genomes/'+gid.value+'/share?fmt='+fmt.value,
   {method:'POST',headers:{Authorization:'Bearer '+token},body:fd});
 if(!r.ok){show(await r.json());return;}
 preview.src=URL.createObjectURL(await r.blob());
 out.textContent='image ready — right-click to save';
}
</script>
"""
