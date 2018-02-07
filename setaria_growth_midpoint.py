# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 10:12:26 2017

@author: Jacob Fredenburg
North Carolina State University
Doherty Lab
"""

import cv2

import argparse

import glob

import math

import csv


'''
    Function to calculate the midpoint between 2 points

    Parameters:
        x1, y1: The coordinates of the first point
        x2, y2: The coordinates of the second point

    Output:
        ymid, xmid: The coordinates of the midpoint
'''
def findMidpoint(x1, y1, x2, y2):
    xmid = (x1 + x2) / 2
    ymid = (y1 + y2) / 2

    return ymid, xmid


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


#parse the commandline arguments to get the path to the folder of images
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--images", required=True, help="full path to input dataset of images")
args = vars(ap.parse_args())

#File name and list of lengths between midpoints
imageInfo = []

#the value used to determine if the RGB values are white
#It is 240 instead of 255 because the edge detection white pixels are not always 255
WHITE_CUTOFF = 240

save = False

ans = input("Save output images? (y/n)\n")

if ans == "y" or ans == "Y" or ans == "yes":
    save = True
elif ans == "n" or ans == "N" or ans == "no":
    save = False
else:
    print("not a valid option. WILL NOT SAVE OUTPUT IMAGES")



for imagePath in glob.glob(args["images"] + "/*.jpg"):

    image = cv2.imread(imagePath)

    print("\nFile: " + imagePath + "\n------------------------------------")

    listOfMidpoints = []

    #start moving through the image
    for y in range(0, image.shape[0]):
        
        #list of white pixels in the row
        whitePixels = []

        lastWhite = False

        #start moving through the row
        for x in range(0, image.shape[1]):

            #Add white pixels to the list
            if image[y][x][0] > WHITE_CUTOFF and image[y][x][1] > WHITE_CUTOFF and image[y][x][2] > WHITE_CUTOFF:
                whitePixels.append([y,x])
                
        #get the number of white pixels in the row
        length = len(whitePixels)
        
        #get the midpoint of the first white pixel and the last white pixel of the row on if there
        #is more than one white pixel in the row
        if length > 1:
            y_mid, x_mid = findMidpoint(whitePixels[length - 1][1], whitePixels[length - 1][0], whitePixels[0][1], whitePixels[0][0])
            listOfMidpoints.append([y_mid,x_mid])
                    


    #Draw the midpoints on the image.
    for point in listOfMidpoints:
        image[int(point[0])][int(point[1])] = [255,255,0]

    leafLength = 0

    #Calculate the distance between the midpoints and add it to length
    for i in range(1, len(listOfMidpoints)):
        point1 = listOfMidpoints[i - 1]
        point2 = listOfMidpoints[i]

        distance = math.sqrt((point2[1] - point1[1])**2 + (point2[0] - point1[0])**2)
        leafLength = leafLength + distance

    print("Length: " + str(leafLength))


    imageInfo.append([getFileName(imagePath), leafLength])

    imagePathSub = imagePath[:-4]

    #Save a new copy of the image with the midpoints drawn
    if save == True:
        print('Saved: ' + str(cv2.imwrite(imagePathSub + "_midpoint.jpg", image, [int(cv2.IMWRITE_JPEG_QUALITY), 100])))


#write the lengths to a output csv file.
try:
    with open(args['images'] + '\\out.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')

        writer.writerow(['File Name'] + ['Length'])
        for l in imageInfo:
            writer.writerow([l[0], l[1]])
except PermissionError:
    print('\nThe output file could not be opened or modified. Check to make sure it is not currently open.')