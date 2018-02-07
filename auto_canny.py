# import the necessary packages
import argparse
import glob
import cv2

    
# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--images", required=True,
	help="full path to input dataset of images")
args = vars(ap.parse_args())

# loop over the images
for imagePath in glob.glob(args["images"] + "/*.jpg"):

	# load the image
	image = cv2.imread(imagePath)
	
	#remove the file extension. Used later to save file
	imagePathSub = imagePath[:-4]
    
    #Apply a median filter to color image
	medianImg = cv2.medianBlur(image, 5)
    
    #Convert median blurred image to grayscale
	gray = cv2.cvtColor(medianImg, cv2.COLOR_BGR2GRAY)
    
    #Gaussian Blur the median blurred image
	blurred = cv2.GaussianBlur(gray, (3, 3), 0)
    
	#Apply Canny edge detection to the Gaussian/Median blurred image
	tightMedian = cv2.Canny(blurred, 75, 250)
	
	#Save image with _medianGaus appended to the file name
	print(cv2.imwrite(imagePathSub + "_medianGaus.jpg", tightMedian, [int(cv2.IMWRITE_JPEG_QUALITY), 100]))
	
