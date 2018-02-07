# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 10:25:24 2017

@author: Jacob Fredenburg
North Carolina State University
Doherty Lab

###############################################################################
###########################List of Constraints#################################
###############################################################################
-Region of interest must be in the middle of the image. Preferred half on one side
half on the other.
-The panicle must be straight up and down
-Only one panicle at a time


###############################################################################
############################Package Versions###################################
###############################################################################
This script was developed using the following versions:

python: 3.6.1 (Anaconda 4.4.0 64-bit)
-numpy: 1.12.1
-cv2: 3.2.0
-sklearn: 0.18.2
-argparse: 1.1


###############################################################################
###################List of Parameters that can be used to tune#################
###############################################################################
-Number of counts that are used to classify a tip or base (TIMES_CUTOFF)
-RMS used to adjust the linear regression (RMS_CUTOFF)

"""
#Provide more powerful arrays
import numpy as np

#Provide image functions
import cv2

#Provide floor and sqrt functions
import math

#Provide functions used to calculate Root Mean Squared Error
from sklearn.metrics import mean_squared_error

#Adds ability to add the need for commandline arguments to the script
import argparse

#Provides functions to read all the contents of a folder
import glob

#Used to time execution, non-essential
import time

#Used to write output file
import csv



#############################################################################
###########################****FUNCTIONS****#################################
#############################################################################

'''
    Function to extract the file name from a path. Used in csv output

    Parameters:
        name: The pathname to the file including the file name

    Output:
        The file name.
'''
def getFileName(name):
    index = -1

    for i in range(len(name) - 1, 0, -1):
        if name[i] == "/" or name[i] == "\\":
            index = i
            break

    return name[index + 1:]

'''
Takes in a start coordinate of a white pixel and "follows" it to the base
add path start

* P P
K W P
K K *

W = white pixel
* = any kind of pixel (white or black)
K = black pixel
P = potential previous pixel

In this case, the W would be classified as a base

This function will follow the white pixels until the scenario above is seen.
When a base is determined it is added to a dictionary keeping track of the number
of times that base has been seen (this is used later to remove extraneous bases).

Parameters:
    startRow: The row of the found white pixel
    startItem: The column of the found white pixel

'''
def followRight(startRow, startItem):
    row = startRow
    item = startItem
    done = False

    while done == False:
            if (row + 1) < image.shape[0] and item < image.shape[1] and item > 1:
                if sum(image[row + 1][item]) < 250 and sum(image[row][item - 1]) < 250 and sum(image[row + 1][item - 1]) < 250:

                    if listOfBases.__contains__([row, item]) == True:

                        countOfPoints[str([row, item])] = countOfPoints[str([row, item])] + 1
                        done = True
                    else:
                        listOfBases.insert(len(listOfBases), [row, item])
                        countOfPoints[str([row, item])] = 1
                        done = True

                elif sum(image[row + 1][item] > 50):
                    row = row + 1
                elif sum(image[row][item - 1]) > 50:
                    item = item - 1
                elif sum(image[row + 1][item - 1] > 50):
                    row = row + 1
                    item = item - 1
                else:
                    return
            else:
                return

'''
The function is the same as the followRight function, but for the left half of the
image.

P P *
P W K
* K K

W = white pixel
* = any kind of pixel (white or black)
K = black pixel
P = potential previous pixel

The base classification scenario is flipped compared to the followRight function.

Parameters:
    startRow: The row of the found white pixel
    startItem: The column of the found white pixel
'''

def followLeft(startRow, startItem):
    row = startRow
    item = startItem
    done = False


    while done == False:
        if (row + 1) < image.shape[0] and item < image.shape[1]:
            if sum(image[row + 1][item]) < 250 and sum(image[row][item + 1]) < 250 and sum(image[row + 1][item + 1]) < 250:

                if listOfBases.__contains__([row, item]) == True:
                    countOfPoints[str([row, item])] = countOfPoints[str([row, item])] + 1
                    done = True
                else:
                    listOfBases.insert(len(listOfBases), [row, item])
                    countOfPoints[str([row, item])] = 1
                    done = True

            elif sum(image[row + 1][item]) > 50:
                   row = row + 1
            elif sum(image[row][item + 1]) > 50:
                item = item + 1
            elif sum(image[row + 1][item + 1]) > 50:
                row = row + 1
                item = item + 1
            else:
                return
        else:
            return
'''
Similar fuction to the follow functions

* K K
P W K
P P *

W = white pixel
* = any kind of pixel (white or black)
K = black pixel
P = potential previous pixel

This function follows white pixels until the above scenario is met. This function
also keeps track of the number of times the tip was encountered.

Parameters:
    startRow: The row of the found white pixel
    startItem: The column of the found white pixel
'''
def tipRight(startRow, startItem):
    row = startRow
    item = startItem
    done = False

    while done == False:
        if row > 0 and row < image.shape[0] - 1 and item < image.shape[1] - 1:
            if sum(image[row - 1][item]) < 250 and sum(image[row][item + 1]) < 250 and sum(image[row - 1][item + 1]) < 250:


                 if listOfTips.__contains__([row, item]) == True:
                     countOfTips[str([row, item])] = countOfTips[str([row, item])] + 1
                     done = True
                 else:
                     listOfTips.insert(len(listOfTips), [row, item])
                     countOfTips[str([row, item])] = 1
                     done = True

            elif sum(image[row - 1][item]) > 50:
                row = row - 1
            elif sum(image[row][item + 1]) > 50:
                item = item + 1
            elif sum(image[row - 1][item + 1]) > 50:
                row = row - 1
                item = item + 1
            else:
                return
        else:
            return


'''
This function has the same purpose as the tipRight function.

K K *
K W P
* P P

W = white pixel
* = any kind of pixel (white or black)
K = black pixel
P = potential previous pixel

The tip classification scenario is flipped compared to the tipRight function.

Parameters:
    startRow: The row of the found white pixel
    startItem: The column of the found white pixel
'''
def tipLeft(startRow, startItem):
    row = startRow
    item = startItem
    done = False

    while done == False:
        if row > 0 and item > 0 and row < image.shape[0] - 1:
            if sum(image[row - 1][item]) < 250 and sum(image[row][item - 1]) < 250 and sum(image[row - 1][item - 1]) < 250:

                if listOfTips.__contains__([row, item]) == True:
                     countOfTips[str([row, item])] = countOfTips[str([row, item])] + 1
                     done = True
                else:
                     listOfTips.insert(len(listOfTips), [row, item])
                     countOfTips[str([row, item])] = 1
                     done = True

            elif sum(image[row - 1][item]) > 50:
                row = row - 1
            elif sum(image[row][item - 1]) > 50:
                item = item -1
            elif sum(image[row - 1][item - 1]) > 50:
                row = row - 1
                item = item - 1
            else:
                return
        else:
            return

'''
This function is responsible for counting the number of pixels between the given
tip and a point given in listOfPoints.

This function will move diagonally before it moves left or right to try and count
the true distance between the tip and base.

Parameters:
    tip: the tip of interest; serves as the start position.
    listOfPoints: a list that contains all the "base" points.

Return:
    pixelCount: The length of the bristle starting a tip
'''
def findBristlesRight(tip, listOfPoints):
    done = False
    row = tip[0]
    item = tip[1]

    pixelCount = 0

    while done == False:
        if (row + 1) < image.shape[0] and item > 0:
            if listOfPoints.__contains__([row, item]) == False:
                if sum(image[row + 1][item - 1]) > 50:

                    row = row + 1
                    item = item - 1
                    pixelCount = pixelCount + 1

                elif sum(image[row][item - 1]) > 50:

                    item = item - 1
                    pixelCount = pixelCount + 1

                elif sum(image[row + 1][item]) > 50:

                    row = row + 1
                    pixelCount = pixelCount + 1

                else:
                    done = True
            else:
                return int(pixelCount)
        else:
            return int(pixelCount)

'''
This function is responsible for counting the number of pixels between the given
tip and a point given in listOfPoints.

This function will move diagonally before it moves left or right to try and count
the true distance between the tip and base.

Parameters:
    tip: the tip of interest; serves as the start position.
    listOfPoints: a list that contains all the "base" points.
'''

def findBristlesLeft(tip, listOfPoints):
    done = False
    row = tip[0]
    item = tip[1]

    pixelCount = 0

    while done == False:
        if (row + 1) < image.shape[0] and item > 0:
            if listOfPoints.__contains__([row, item]) == False:
                if sum(image[row + 1][item + 1]) > 50:

                    row = row + 1
                    item = item + 1
                    pixelCount = pixelCount + 1

                elif sum(image[row][item + 1]) > 50:

                    item = item + 1
                    pixelCount = pixelCount + 1

                elif sum(image[row + 1][item]) > 50:

                    row = row + 1
                    pixelCount = pixelCount + 1

                else:
                    done = True
            else:
                return int(pixelCount)
        else:
            return int(pixelCount)

'''
This functions sorts the length of the bristles and returns parallel lists one
containing tip coordinates and the other containing the lengths.


Parameters:
    lengths: a dictionary of the string version of the tip coordinates as the
    key and the length as the value.

Return:
    sorted_Bristle_Tips: A list containing the sorted tip coordinates.
    sorted_Bristle_Lengths: A list containing the sorted lengths of the bristles


The bristle starting at sorted_Bristle_Tips[0] has length sorted_Bristle_Lengths[0]

'''
def determineLongestBristles(lengths):

    sorted_Bristle_Tips = []
    sorted_Bristle_Lengths = []

    for k,v in lengths.items():
        if v != None:
            sorted_Bristle_Tips.insert(0,k)
            sorted_Bristle_Lengths.insert(0,v)


    for i in range(0, len(sorted_Bristle_Lengths)):
        for j in range(0, len(sorted_Bristle_Lengths)):
            if sorted_Bristle_Lengths[i] > sorted_Bristle_Lengths[j]:
                temp = sorted_Bristle_Lengths[i]
                sorted_Bristle_Lengths[i] = sorted_Bristle_Lengths[j]
                sorted_Bristle_Lengths[j] = temp

                temp = sorted_Bristle_Tips[i]
                sorted_Bristle_Tips[i] = sorted_Bristle_Tips[j]
                sorted_Bristle_Tips[j] = temp

    return sorted_Bristle_Tips, sorted_Bristle_Lengths


#############################################################################
#########################****START BODY CODE****#############################
#############################################################################


ap = argparse.ArgumentParser()
ap.add_argument("-i", "--images", required=True, help="full path to input dataset of images")
args = vars(ap.parse_args())

save = False

ans = input("Save output images? (y/n)\n")

if ans == "y" or ans == "Y" or ans == "yes":
    save = True
elif ans == "n" or ans == "N" or ans == "no":
    save = False
else:
    print("not a valid option. WILL NOT SAVE OUTPUT IMAGES")


t0 = time.time()

#list used to hold all the features of the each panicle in each image
features = []

#Process each image in the destination folder
for imagePath in glob.glob(args["images"] + "/*.jpg"):

    #List used to hold the total height (y-value), bristle density, and the max width
    #to be written out to a csv file
    listOfFeatures = []

    print("\nFile: " + imagePath + "\n------------------------------------")
    listOfFeatures.append(getFileName(imagePath))

    #Load image
    #Image is stored as a 3d list in this format: (Y, X, BGR vals)
    image = cv2.imread(imagePath)

    #Variables used to iterate over the pixels in the image
    rowCount = 0
    itemCount = 0

    #List to keep track of the postitions of the bases and tips encountered
    listOfBases = []
    listOfTips = []


    #Dictionaries to keep track of how many times each base or tip in encountered
        #key: coordinates of the base/tip
        #value: number of times the base/tip was encountered
    countOfPoints = {}
    countOfTips = {}

    #Calculation of the middle of the image. Gives 'right' side of image the fractional
    #pixel if there is one
    middle = math.floor(image.shape[1] / 2)

    #%%
    #############################################################################
    ###########################****FIND BASES****################################
    #############################################################################

    #Iterate through the image and if a white pixel is found call 'follow' to
    #find tips and bases
    firstWhite = 0
    foundFirst = False
    for row in image:
        itemCount = 0
        for item in row:
            if item[0] >= 250 and item[1] >= 250 and item[2] >= 250 and itemCount >= middle:
                if foundFirst == False:
                    firstWhite = rowCount
                    foundFirst = True

                followRight(rowCount, itemCount)
                tipRight(rowCount, itemCount)
            elif item[0] >= 250 and item[1] >= 250 and item[2] >= 250 and itemCount < middle:
                if foundFirst == False:
                    firstWhite = rowCount
                    foundFirst = True

                followLeft(rowCount, itemCount)
                tipLeft(rowCount, itemCount)

            itemCount = itemCount + 1
        rowCount = rowCount + 1


    #Change the color of the pixel for each base and tip
    for point in listOfBases:
        #Blue
        image.itemset((point[0], point[1], 0), 0)

        #Green
        image.itemset((point[0], point[1], 1), 255)

        #Red
        image.itemset((point[0], point[1], 2), 127)

    for tip in listOfTips:
        #Blue
        image.itemset((tip[0], tip[1], 0), 220)

        #Green
        image.itemset((tip[0], tip[1], 1), 66)

        #Red
        image.itemset((tip[0], tip[1], 2), 244)



    #Remove bases and tips that were not encountered a specific set of times
    TIMES_CUTOFF = 20

    for k, v in countOfPoints.items():
        if v < TIMES_CUTOFF:
            k = k[:-1]
            k = k[1:]
            x, y = k.split(", ")
            listOfBases.remove([int(x), int(y)])
            image[int(x)][int(y)] = [255, 255, 255]

    for k, v in countOfTips.items():
        if v < TIMES_CUTOFF:
            k = k[:-1]
            k = k[1:]
            x, y = k.split(", ")
            listOfTips.remove([int(x), int(y)])
            image[int(x)][int(y)] = [255, 255, 255]




    #%%
    #############################################################################
    #####################****SEPARATE POINTS FOR LINE REG****####################
    #############################################################################

    #Lists to keep track of bases and tips from each side of the image separately
    listOfBasesRight = []
    listOfBasesLeft = []
    listOfTipsRight = []
    listOfTipsLeft = []

    #Separate the right side of the image
    for point in listOfBases:
        if point[1] >= middle:
            listOfBasesRight.append(point)
        elif point[1] < middle:
            listOfBasesLeft.append(point)
        else:
            print(str(point) + "was not placed in a linear reg category")

    #Separate the left side of the image
    for tip in listOfTips:
        if tip[1] >= middle:
            listOfTipsRight.append(tip)
        elif tip[1] < middle:
            listOfTipsLeft.append(tip)
        else:
            print(str(tip) + " was not placed in a linear reg category")

    #Dictionary to to hold all the lengths of the bristles
    bristleLengths = {}

    #Add the right side of the image to the dictionary
    for tip in listOfTipsRight:
        bristleLengths[str(tip)] = findBristlesRight(tip, listOfBases)

    #Add the left side of the image to the dictionary
    for tip in listOfTipsLeft:
        bristleLengths[str(tip)] = findBristlesLeft(tip, listOfBases)


    #Get the sorted parallel lists for longest bristles
    longest_Tips, longest_Lengths = determineLongestBristles(bristleLengths)



    #Separate point points into x and y for plotting and lin reg
    listOfBasesRightX = []
    listOfBasesRightY = []

    listOfBasesLeftX = []
    listOfBasesLeftY = []

    listOfTipsRightX = []
    listOfTipsRightY = []

    listOfTipsLeftX = []
    listOfTipsLeftY = []

    for point in listOfBasesRight:
        listOfBasesRightX.append(point[1])
        listOfBasesRightY.append(point[0])

    for point in listOfBasesLeft:
        listOfBasesLeftX.append(point[1])
        listOfBasesLeftY.append(point[0])

    for tip in listOfTipsRight:
        listOfTipsRightX.append(tip[1])
        listOfTipsRightY.append(tip[0])

    for tip in listOfTipsLeft:
        listOfTipsLeftX.append(tip[1])
        listOfTipsLeftY.append(tip[0])

    #%%
    #############################################################################
    ##########################****LINEAR REGRESSION****##########################
    #############################################################################

    #Rotation of the image is necessary to the linear regression to work
    #Redfining the separate x and y lists to make code easier to understand

    listOfBasesRightX_ROTATED = listOfBasesRightY
    listOfBasesRightY_ROTATED = listOfBasesRightX

    listOfBasesLeftX_ROTATED = listOfBasesLeftY
    listOfBasesLeftY_ROTATED = listOfBasesLeftX

    listOfTipsRightX_ROTATED = listOfTipsRightY
    listOfTipsRightY_ROTATED = listOfTipsRightX

    listOfTipsLeftX_ROTATED = listOfTipsLeftY
    listOfTipsLeftY_ROTATED = listOfTipsLeftX


    #Get quadratic regression coefficients
    Line_Reg_Right = np.polyfit(listOfBasesRightX_ROTATED, listOfBasesRightY_ROTATED, deg = 2)
    Line_Reg_Left = np.polyfit(listOfBasesLeftX_ROTATED, listOfBasesLeftY_ROTATED, deg = 2)


    #Cutoff for the RMS to compare to. If the calculated RMS is greater than this
    #the bases is discarded.
    RMS_CUTOFF = 100

    #############################################################################
    #Adjust Left Linear Regression
    #############################################################################
    index = 0
    pointCount = 0
    removed = 0
    for p in listOfBasesLeftX_ROTATED:

        #Gets actual y(rotated) value
        y = listOfBasesLeftY_ROTATED[index]

        #Gets the predicted y(rotated) value
        pred = Line_Reg_Left[2] + Line_Reg_Left[1] * p + Line_Reg_Left[0] * (p**2)

        #calculates the RSME
        rms = math.sqrt(mean_squared_error([y], [pred]))

        if rms > RMS_CUTOFF:
            listOfBasesLeftX_ROTATED[index] = -1
            listOfBasesLeftY_ROTATED[index] = -1


        index = index + 1
        pointCount = pointCount + 1

    done = False
    while done == False:
        try:
            listOfBasesLeftX_ROTATED.remove(-1)
            listOfBasesLeftY_ROTATED.remove(-1)
        except ValueError:
            done = True


    #############################################################################
    #Adjust Right Linear Regression
    #############################################################################
    index = 0
    pointCount = 0
    removed = 0
    for p in listOfBasesRightX_ROTATED:
        #Gets actual y(rotated) value
        y = listOfBasesRightY_ROTATED[index]

        #Gets the predicted y(rotated) value
        pred = Line_Reg_Right[2] + Line_Reg_Right[1] * p + Line_Reg_Right[0] * (p**2)

        #calculates the RSME
        rms = math.sqrt(mean_squared_error([y], [pred]))

        if rms > RMS_CUTOFF:
            listOfBasesRightX_ROTATED[index] = -1
            listOfBasesRightY_ROTATED[index] = -1

        index = index + 1
        pointCount = pointCount + 1

    done = False
    while done == False:
        try:
            listOfBasesRightX_ROTATED.remove(-1)
            listOfBasesRightY_ROTATED.remove(-1)
        except ValueError:
            done = True

    #Adjust the Lin Reg line
    Line_Reg_Right_Ad = np.polyfit(listOfBasesRightX_ROTATED, listOfBasesRightY_ROTATED, deg = 2)
    Line_Reg_Left_Ad = np.polyfit(listOfBasesLeftX_ROTATED, listOfBasesLeftY_ROTATED, deg = 2)



    #%%
    #############################################################################
    ##################****CALCULATE PANICLE HEIGHT****###########################
    #############################################################################

    #Lists used to hold the x and y values of the intersection points between
    #the regression lines
    x_Values = []
    y_Values = []

    #Calculate the a, b, and c term for the quadratic equation
    a = (Line_Reg_Left_Ad[0] - Line_Reg_Right_Ad[0])
    b = (Line_Reg_Left_Ad[1] - Line_Reg_Right_Ad[1])
    c = (Line_Reg_Left_Ad[2] - Line_Reg_Right_Ad[2])

    #Calculate the determinant of the quadratic equation
    det = (b**2) - 4*(a*c)

    #As long as the det is not negative, solve for the intersection points
    #otherwise use the image dimensions as the x-values. The x-values are used
    #later
    if det > 0:

        x1 = (-b - math.sqrt( (b**2)-(4*(a*c)))) / (2*a)
        x2 = (-b + math.sqrt((b**2)-(4*(a*c))))/(2*a)

        x_Values.append(math.floor(x1))
        x_Values.append(math.floor(x2))



        y1 = Line_Reg_Left_Ad[2] + Line_Reg_Left_Ad[1] * x_Values[0] + Line_Reg_Left_Ad[0] * (x_Values[0]**2)
        y2 = Line_Reg_Left_Ad[2] + Line_Reg_Left_Ad[1] * x2 + Line_Reg_Left_Ad[0] * (x2**2)

        y_Values.append(math.floor(y1))
        y_Values.append(math.floor(y2))

        print("Total Height: " + str(x_Values[0]))
        listOfFeatures.append(x_Values[0])

        distance = math.sqrt((x_Values[1] - x_Values[0])**2 + (y_Values[1] - y_Values[0])**2)

        print("Panicle Length: " + str(distance))
        listOfFeatures.append(distance)

        #Use the bottom of the image and the lesser intersection point for the total height of the plant
        bristleDensity = len(longest_Tips) / (image.shape[0] - x_Values[0])
        print("Bristle Density: " + str(bristleDensity))

        listOfFeatures.append(bristleDensity)

    else:
        print("The regression lines never intersect. Using the row of the first white pixel instead")

        print("Total Height: " + str(firstWhite))
        listOfFeatures.append(firstWhite)

        print("Panicle Length: " + str(firstWhite))
        listOfFeatures.append(firstWhite)

         #Use the firstWhite pixel and the bottom of the image as total height of the plant
        bristleDensity = len(longest_Tips) / (image.shape[0] - firstWhite)

        print("Bristle Density: " + str(bristleDensity))
        listOfFeatures.append(bristleDensity)

        x_Values = [firstWhite, image.shape[0]]


    #%%
    #############################################################################
    ###################****FIND MAX PANICLE WIDTH****############################
    #############################################################################

    maxWidth = 0
    atX = 0
    atY = [0, 0]

    #Iterate between the x-values calculated above to determine the max width
    for val in range(x_Values[0], x_Values[1]):

        #Calculate the y values for both intersection lines.
        rightVal = Line_Reg_Right_Ad[2] + Line_Reg_Right_Ad[1] * val + Line_Reg_Right_Ad[0] * (val**2)
        leftVal = Line_Reg_Left_Ad[2] + Line_Reg_Left_Ad[1] * val + Line_Reg_Left_Ad[0] * (val**2)

        if abs(rightVal - leftVal) > maxWidth:
            maxWidth = abs(rightVal - leftVal)
            atX = val
            atY[0] = rightVal
            atY[1] = leftVal

    print("Max width: " + str(maxWidth))
    listOfFeatures.append(maxWidth)
    x1 = [atX, atX]



    #%%
    #############################################################################
    ###############################****PLOTTING****##############################
    #############################################################################

    x = np.linspace(0, image.shape[0], image.shape[0])
    x_round = []
    x_round.extend(range(0, image.shape[0]))


    ###############################################################################
    # PLOT THE REGRESSION LINES ON THE FOLLOWED IMAGE
    ###############################################################################
    y = Line_Reg_Left_Ad[2] + Line_Reg_Left_Ad[1] * x + Line_Reg_Left_Ad[0] * (x**2)
    y_round = np.round(y)


    ignoredPixels = []
    for i in range(0, image.shape[0]):
        try:
            #Blue
            image.itemset((x_round[i], int(y_round[i]), 0), 0)

            #Green
            image.itemset((x_round[i], int(y_round[i]), 1), 255)

            #Red
            image.itemset((x_round[i], int(y_round[i]), 2), 255)
        except IndexError:

            #Ignore the pixel
            ignoredPixels.append([x_round[i], y_round[i]])


    y = Line_Reg_Right_Ad[2] + Line_Reg_Right_Ad[1] * x + Line_Reg_Right_Ad[0] * (x**2)
    y_round = np.round(y)


    for i in range(0, image.shape[0]):
        try:
            #Blue
            image.itemset((x_round[i], int(y_round[i]), 0), 255)

            #Green
            image.itemset((x_round[i], int(y_round[i]), 1), 255)

            #Red
            image.itemset((x_round[i], int(y_round[i]), 2), 0)
        except IndexError:

            #Ignore the pixel
            ignoredPixels.append([x_round[i], y_round[i]])


    features.append(listOfFeatures)

    imagePathSub = imagePath[:-4]

    #Save the "followed" image
    if save == True:
        print( "Saved: " + str(cv2.imwrite(imagePathSub + "_followed.jpg", image, [int(cv2.IMWRITE_JPEG_QUALITY), 100])))


#Opens a csv file named out in the location of the input images

#Writes the total height, bristle density, max width to the file for easier data
#entry
try:
    with open(args['images'] + '\\out.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['File'] + ['Total Height'] + ['Panicle Length'] + ['Bristle Density'] + ['Max Width'])
        for l in features:
            writer.writerow([l[0], l[1], l[2], l[3], l[4]])
except PermissionError:
    print('\nThe output file could not be opened or modified. Check to make sure it is not currently open.')


t1 = time.time()

print("It took " + str(t1 - t0) + " seconds")
