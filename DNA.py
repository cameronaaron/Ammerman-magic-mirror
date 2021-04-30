#!/usr/bin/env python
# coding: utf-8

# #!/usr/bin/env python3

# In[1]:


#get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
import sys
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
from tkinter import *
from tkinter.ttk import *
from tkinter.filedialog import askopenfile
import tkinter as tk
from tkinter import filedialog
Fileexistprompt=input('Do you have a 23andme file?').lower()
if Fileexistprompt == 'yes':
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename()
    data = pd.read_csv(file_path, sep='\t', dtype={'rsid':'str', 'chromosome':'object', 'position':'int', 'genotype':'str'}, comment='#')
    btn = Button(root, text ='Open', command = lambda:open_file())
    btn.pack(side = TOP, pady = 10)
elif Fileexistprompt == 'no':
    data = pd.read_csv("dnafiles/genomee.txt", sep='\t', dtype={'rsid':'str', 'chromosome':'object', 'position':'int', 'genotype':'str'}, comment='#')


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
text = genes_to_display.summary.values
stop_words = ["Abscess",
"Acute Radiation Sickness",
"Alzheimers",
"Alzheimer",
"likely",
"met",
"disease",
"Anthrax",
"Appendicitis",
"Allergy",
"Arthritis",
"Aseptic meningitis",
"Asthma",
"Astigmatism",
"Atherosclerosis",
"B",
"Bacterial meningitis",
"Beriberi",
"Black Death",
"Botulism",
"Breast cancer",
"Bronchitis",
"Brucellosis",
"Bubonic plague",
"Bunion",
"Boil",
"C",
"Campylobacter infection",
"Cancer",
"Candidiasis",
"Carbon monoxide poisoning",
"Coeliac disease",
"Cerebral palsy",
"Chagas disease",
"Chickenpox",
"Chlamydia",
"Chlamydia trachomatis",
"Cholera",
"Chordoma",
"Chorea",
"Chronic fatigue syndrome",
"Circadian rhythm sleep disorder",
"Colitis",
"Common cold",
"Condyloma",
"Congestive heart disease",
"Coronary heart disease",
"COVID-19",
"Cowpox",
"Crohn's Disease",
"Coronavirus",
"D",
"Dengue Fever",
"Diabetes mellitus",
"Diphtheria",
"Dehydration",
"Dysentery",
"E",
"Ear infection",
"Ebola",
"Encephalitis",
"Emphysema",
"Epilepsy",
"Erectile dysfunction",
"F",
"Fibromyalgia",
"Foodborne illness",
"G",
"normal",
"prp",
"cjd",
"Normal",
"val",
"vcjd",
"Gangrene",
"Gastroenteritis",
"Genital herpes",
"GERD",
"Goitre",
"Gonorrhea",
"H",
"Heart disease",
"Hepatitis A",
"Hepatitis B",
"Hepatitis C",
"Hepatitis D",
"Hepatitis E",
"Histiocytosis (childhood cancer)",
"HIV",
"Human papillomavirus",
"Huntington's disease",
"Hypermetropia",
"Hyperopia",
"Hyperthyroidism",
"Hypothyroid",
"Hypotonia",
"I",
"Impetigo",
"Infertility",
"Influenza",
"Interstitial cystitis",
"Iritis",
"Iron-deficiency anemia",
"Irritable bowel syndrome",
"Ignious Syndrome",
"Intestine ache",
"Intestine Gas",
"Intestine disease",
"Upset Intestine",
"J",
"Jaundice",
"K",
"Keloids",
"Kuru",
"Kwashiorkor",
"Kidney stone disease",
"L",
"Laryngitis",
"Lead poisoning",
"Legionellosis",
"Leishmaniasis",
"Leprosy",
"Leptospirosis",
"Listeriosis",
"Leukemia",
"Lice",
"Loiasis",
"Lung cancer",
"Lupus erythematosus",
"Lyme disease",
"Lymphogranuloma venereum",
"Lymphoma",
"Limbtoosa",
"M",
"Mad cow disease",
"Malaria",
"Marburg fever",
"Measles",
"Melanoma",
"Metastatic cancer",
"Meniere's disease",
"Meningitis",
"Migraine",
"Mononucleosis",
"Multiple myeloma",
"Multiple sclerosis",
"Mumps",
"Muscular dystrophy",
"Myasthenia gravis",
"Myelitis",
"Cancer",
"Cancer'",
"Increased",
"'",
"Myoclonus",
"Myopia",
"Myxedema",
"Morquio Syndrome",
"Mattticular syndrome",
"Mononucleosis",
"N",
"Neoplasm",
"Non-gonococcal urethritis",
"Necrotizing Fasciitis",
"Night blindness",
"O",
"Obesity",
"Osteoarthritis",
"Osteoporosis",
"Otitis",
"P",
"Cancer's",
"Palindromic rheumatism",
"Paratyphoid fever",
"Parkinson's disease",
"Pelvic inflammatory disease",
"Peritonitis",
"Periodontal disease",
"Pertussis",
"Phenylketonuria",
"Plague",
"Poliomyelitis",
"Porphyria",
"Progeria",
"Prostatitis",
"Psittacosis",
"Psoriasis",
"Pubic lice",
"Pulmonary embolism",
"Pilia",
"pneumonia",
"Q",
"Q fever",
"Ques fever",
"R",
"Rabies",
"Repetitive strain injury",
"Rheumatic fever",
"Rheumatic heart",
"Rheumatism",
"Rheumatoid arthritis",
"Rickets",
"Rift Valley fever",
"Rocky Mountain spotted fever",
"Rubella",
"S",
"Salmonellosis",
"Scabies",
"Scarlet fever",
"Sciatica",
"Scleroderma",
"Scrapie",
"Scurvy",
"Sepsis",
"Septicemia",
"SARS",
"Shigellosis",
"Shin splints",
"Shingles",
"Sickle-cell anemia",
"Siderosis",
"SIDS",
"Silicosis",
"Smallpox",
"Stevens–Johnson syndrome",
"Stomach flu",
"Stomach ulcers",
"Strabismus",
"Strep throat",
"Streptococcal infection",
"Synovitis",
"Syphilis",
"Swine influenza",
"Schizophrenia",
"Stomach Gas",
"Stomach Ache",
"stomach Disease",
"Kids Stomach Ache",
"Upset Stomach",
"T",
"Taeniasis",
"Tay-Sachs disease",
"Tennis elbow",
"Teratoma",
"Tetanus",
"Thalassaemia",
"Thrush",
"Thymoma",
"Tinnitus",
"Tonsillitis",
"Tooth decay",
"Toxic shock syndrome",
"Trichinosis",
"Trichomoniasis",
"Trisomy",
"Tuberculosis",
"Tularemia",
"Tungiasis",
"Typhoid fever",
"Typhus",
"Tumor",
"U",
"Ulcerative colitis",
"Ulcers",
"Uremia",
"Urticaria",
"Uveitis",
"UTI'S",
"V",
"risk",
"Varicella",
"Varicose veins",
"Vasovagal syncope",
"Vitiligo",
"Von Hippel-Lindau disease",
"Viral fever",
"Viral meningitis",
"W",
"Warkany syndrome",
"Warts",
"Watkins",
"Y",
"Yellow fever",
"Yersiniosis"] + list(STOPWORDS)

wordcloud = WordCloud(
    width = 800,
    height = 600,
    mode = 'RGBA',
    background_color = None,
    stopwords = stop_words).generate(str(text))

fig = plt.figure(
    figsize = (90, 80),
    edgecolor = 'k')
plt.imshow(wordcloud, interpolation = 'bilinear')
plt.axis('off')
plt.tight_layout(pad=0)
plt.show()
wordcloud.to_file("wordcloud.png")
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
