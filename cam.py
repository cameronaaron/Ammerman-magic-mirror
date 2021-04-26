import cv2

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
