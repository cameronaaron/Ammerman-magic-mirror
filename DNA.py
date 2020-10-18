#!/usr/bin/env python
# coding: utf-8

# #!/usr/bin/env python3

# In[1]:


#get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
sns.set_style('darkgrid')
sns.color_palette('Spectral')
import matplotlib.pyplot as plt
import numpy as np
import requests
import pandas as pd
import re
import matplotlib.pyplot as plt
import cv2
from wordcloud import WordCloud, STOPWORDS
from PIL import Image
from lineage import Lineage


# In[2]:


data = pd.read_csv('genomee.txt', sep='\t', dtype={'rsid':'str', 'chromosome':'object', 'position':'int', 'genotype':'str'}, comment='#')


# In[3]:


print(data)


# In[4]:


df = pd.DataFrame(data)


# In[5]:


df.head(25)


# In[6]:


df.nunique()


# In[7]:

'''
duplicates = df[df.duplicated(subset='position')]
display(duplicates.head())
display(duplicates.info())
'''

# In[8]:


df = df[df.chromosome != 'Y']
df.info()


# In[9]:


df['chromosome'].unique()


# In[10]:


df['chromosome'] = df['chromosome'].apply(lambda x: re.sub(r'X', r'23', x))
df['chromosome'] = df['chromosome'].apply(lambda x: re.sub(r'MT', r'24', x))


# In[11]:


df['chromosome'] = df['chromosome'].apply(lambda x: int(x))


# In[12]:


chromosome_dict = {1:'1', 2:'2', 3:'3', 4:'4', 5:'5', 6:'6', 7:'7', 8:'8', 9:'9', 10:'10', 11:'11', 12:'12', 13:'13',
                  14:'14', 15:'15', 16:'16', 17:'17', 18:'18', 19:'19', 20:'20', 21:'21', 22:'22', 23:'X', 24:'MT'}


# In[13]:


print(chromosome_dict)
df.info()


# In[14]:


genotype_na = df[df.genotype == '--']
len(genotype_na)


# In[15]:


df[df.chromosome == 1].info()


# In[16]:


df.rename({' rsid': 'rsid'}, axis='columns', inplace=True)


# In[17]:


rsid_per_chromosome_series = df.groupby('chromosome')['rsid'].count()
rsid_per_chromosome_series.columns = ['chromosome', 'count']


# In[18]:


rsid_per_chromosome_series.plot.barh(figsize=(16,9), fontsize=15)
plt.show()
plt.savefig("counts.png")


# In[19]:


snp_df = pd.read_csv('result.csv')
snp_df.head()


# In[20]:


snp_df['genotype'] = snp_df['Unnamed: 0'].apply(lambda x: re.sub(r'.*([AGCT]);([AGCT])\)', r'\1\2', x))


# In[21]:


snp_df.head()


# In[22]:


new_cols = ['rsid', 'magnitude', 'summary', 'genotype']
snp_df.columns = new_cols


# In[23]:


snp_df['rsid'] = snp_df['rsid'].map(lambda x : x.lower())
snp_df['rsid'] = snp_df['rsid'].map(lambda x : re.sub(r'([a-z]{1,}[\d]+)\([agct];[agct]\)', r'\1', x))


# In[24]:


snp_df.head()


# In[25]:


snp_df.info()


# In[26]:


snp_df.isna().any()


# In[27]:


new_df = snp_df.merge(df, how='inner', on=['rsid', 'genotype'], suffixes=('_SNPedia', '_myDNA'))


# In[28]:


new_df.head(1000000)


# In[29]:


genes_to_display = new_df[new_df.magnitude > 2]


# In[30]:


genes_to_display


# In[31]:


print (genes_to_display.summary.values)


# In[32]:


from wordcloud import WordCloud, STOPWORDS
import matplotlib.pyplot as plt
text = genes_to_display.summary.values
wordcloud = WordCloud(
    width = 800,
    height = 600,
    mode = 'RGBA',
    background_color = None,
    stopwords = STOPWORDS).generate(str(text))
fig = plt.figure(
    figsize = (90, 80),
    edgecolor = 'k')
plt.imshow(wordcloud, interpolation = 'bilinear')
plt.axis('off')
plt.tight_layout(pad=0)
plt.show()
wordcloud.to_file("wordcloud.png")

# In[33]:

'''
from PIL import Image

img = Image.open("wordcloud.png")
img = img.convert("RGBA")
datas = img.getdata()

newData = []
for item in datas:
    if item[0] == 255 and item[1] == 255 and item[2] == 255:
        newData.append((255, 255, 255, 0))
    else:
        if item[0] > 150:
            newData.append((0, 0, 0, 255))
        else:
            newData.append(item)
            print(item)


img.putdata(newData)
img.save("wordcloud.png", "PNG")
'''

# In[ ]:


deviceId = 0
faceCascade = cv2.CascadeClassifier("cascadeFiles/haarcascade_frontalface_default.xml")
noseCascade = cv2.CascadeClassifier("reference/haarcascade_mcs_nose.xml")
eyeCascade = cv2.CascadeClassifier("cascadeFiles/haarcascade_eye.xml")
leftEyeCascade = cv2.CascadeClassifier("cascadeFiles/haarcascade_lefteye_2splits.xml")


imgHat = cv2.imread('wordcloud.png',-1)

orig_mask = imgHat[:,:,3]

orig_mask_inv = cv2.bitwise_not(orig_mask)

imgHat = imgHat[:,:,0:3]
origHatHeight, origHatWidth = imgHat.shape[:2]

cv2.namedWindow("Live Feed", 0)
cv2.setWindowProperty("Live Feed", 0, 1)

video_capture = cv2.VideoCapture(deviceId)

while(cv2.waitKey(30) != 27):
    ret, frame = video_capture.read()
    height,width,_ = frame.shape
    overlayed = frame
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

    faces = faceCascade.detectMultiScale(
        gray,
        scaleFactor=1.1,
        minNeighbors=5,
        minSize=(30, 30),
        flags=cv2.CASCADE_SCALE_IMAGE
    )

    x = 0
    y = 0
    w = 0
    h = 0

    for (tx, ty, tw, th) in faces:
        if tw*th > w*h:
            x = tx
            y = ty
            w = tw
            h = th

    #print x
    #print y
    #print w
    #print h
    hatWidth = 0
    hatHeight = 0
    if (w != 0) and (h != 0):
        face = cv2.rectangle(frame,(x,y),(x+w,y+h),(255,0,0),2)
        hatWidth = (int)(w * 3 / 2)
        hatHeight = (int)(hatWidth * origHatHeight / origHatWidth)

        x1 = (int)(x - (hatWidth/4))
        x2 = (int)(x + w + (hatWidth/4))
        y1 = (int)(y - (hatHeight*3/4))
        y2 = (int)(y + h+ (hatHeight/4))

        if x1 < 0:
            x1 = 0
        if x2 > width:
            x2 = width
        if y1 < 0:
            y1 = 0
        if y2 > height:
            y2 = height

        hatHeight = (int)(y2 - y1)
        hatWidth = (int)(x2 - x1)

        hat = cv2.resize(imgHat, (hatWidth,hatHeight), interpolation = cv2.INTER_AREA)
        mask = cv2.resize(orig_mask, (hatWidth,hatHeight), interpolation = cv2.INTER_AREA)
        mask_inv = cv2.resize(orig_mask_inv, (hatWidth,hatHeight), interpolation = cv2.INTER_AREA)

        roi = frame[y1:y2, x1:x2]
        try:
            roi_bg = cv2.bitwise_and(roi, roi, mask=mask_inv)
            roi_fg = cv2.bitwise_and(hat, hat, mask=mask)
            dst = cv2.add(roi_bg,roi_fg)
            frame[y1:y2, x1:x2] = dst


            roi_gray_m = gray[y:y+h, x:x+w]
            roi_color_m = frame[y:y+h, x:x+w]


        finally:
            cropy1 = y + (h/2) - ((x2-x1)*(float(2)/3))
            cropy2 = y + (h/2) + ((x2-x1)*(float(2)/3))

            if cropy1 < 0:
                cropy1 = 0
            if(cropy1 >= cropy2):
                cropy1 = cropy2-1
            if cropy2 >height:
                cropy2=height
            print ("cropy2: ")
            print (cropy2)
            print ("cropy1: ")
            print (cropy1)
            print ("x1: ")
            print (x1)
            print ("x2: ")
            print (x2)
            #small = cv2.resize(frame[5:100, 5:x2], (hatWidth,hatHeight), fx=0.5, fy=0.5)
            overlayed = frame[int(cropy1):int(cropy2), x1:x2]


    cv2.imshow("Live Feed", overlayed)

    if cv2.waitKey(1) & 0xFF == ord('q'):
        break
video_capture.release()
cv2.destroyAllWindows()

'''
# In[ ]:


l = Lineage()


# In[ ]:


Mom = l.create_individual('mom', 'mom.txt')


# In[ ]:


Mom


# In[ ]:


Mom.build


# In[ ]:


Me = l.create_individual('me', 'genomeeline.txt')


# In[ ]:


Me.build


# In[ ]:


discordant_snps = l.find_discordant_snps(Me, Mom, save_output=True)


# In[ ]:


len(discordant_snps.loc[discordant_snps['chrom'] != 'MT'])


# In[ ]:


results = l.find_shared_dna([Mom, Me], cM_threshold=0.75, snp_threshold=1100)


# In[ ]:


sorted(results.keys())


# In[ ]:


len(results['one_chrom_shared_dna'])


# In[ ]:


results1 = l.find_shared_dna([Mom, Me], shared_genes=True)


# In[ ]:


len(results1['two_chrom_shared_genes'])


# In[ ]:


from arv import load, unphased_match as match

genome = load("genomee.txt")

print("You are {gender}. You are {athletic}. you tend to {sneezesun} when it is too sunny outside. you have {color} eyes,  {hair} hair. You are likely {bodytype} and still struggles with {OCD}. you {likelyhiv}. you are {social} and {popularity}, and {organization}. You are an extremely {drive}. You are quite a {empathy}. You are {likelylynch}. You are {lactoseint}."
.format(
  gender     = "man" if genome.y_chromosome else "woman",
  athletic   = "Incredibly athletic and likely a sprinter" if genome["rs1815739"] == "CC" else "not athletic",
  sneezesun  ="Sneeze" if genome["rs10427255"] == "CC" else "not sneeze",
  social     ="Very social " if genome["rs53576"] == "AA" or "AG" else "not very social",
  hair       ="curly" if genome["rs17646946"] == "GG" else "not curly",
  bodytype   ="Muscular" if genome["rs1815739"] == "CC" else "not Muscular",
  organization ="organized" if genome["rs25532"] == "CC" or "CT" else "Not organizations",
  likelyhiv  ="Hiv resistant" if genome["i3003626"] == "DD" else "are not hiv resistant ",
  empathy     ="Very empathetic " if genome["rs53576"] == "AA" or "AG" else "not very empathetic",
  popularity ="very popular" if genome["rs53576"] == "AA" or "AG" else "not very popular",
  drive      ="very driven" if genome["rs1815739"] == "AA" else "not very driven",
  OCD        ="OCD" if genome["rs25532"] == "CC" or "CT" else "Not ocd",
  likelylynch="Have Lynch syndrome" if genome["rs63750875"] == "CC" or "CG" else "dont have lynch syndrome",
  lactoseint ="lactose intolerent" if genome["rs4988235"] == "CC"  else "not lactose intolerant",
  complexion = "light" if genome["rs1426654"] == "AA" else "dark",
  color      = match(genome["rs12913832"], {"AA": "brown",
                                            "AG": "brown or green",
                                            "GG": "blue"})))
'''
