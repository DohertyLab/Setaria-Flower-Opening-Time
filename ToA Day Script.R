
##how to calculate daytime image flower opening time
# InDir is a directory of all the panicle images for one day
# By is the sequence of images to take from the directory with images (ie By=30 would take every 30th image in the directory)
# mx.model is the MXnet model used for classification
# iter.mx is the iteriation of the MXnet model used
# resize_factor is the numeric between 0-1 that reduces the size of the image by the factor.
Day_spikes=ReadMulti_DayPanicle(InDir = TheDir,By = 30,
                                mx.model = '/mnt/scratch/jigar/flower_opening_learning/Inception_Small-open_closed', 
                                iter.mx = 85, resize_factor = .75)


#Function for daytime images Flower openning time
ReadMulti_DayPanicle<-function(InDir, By=1, first_jpg, mx.model=NULL, iter.mx=NULL,  resize_factor=.7)
{
  
  #load images into stack
  paste("Select Panicle Area in Plot Area")  
  AllFiles=dir(path = InDir, pattern = ".jpg$", full.names = T)
  AllFiles<-sort(AllFiles)
  AllFiles<-AllFiles[seq(1,length(AllFiles),by = By)]
  names(AllFiles)<-AllFiles
  
  #Adjust the images by croping out the panicle from the image.
  timg=as.raster(load.image(file = AllFiles[1]))
  plot(timg)
  bg.crop = base::as.vector(extent(select(raster(x = AllFiles[1], band = 2))))
  timg=timg[(dim(timg)[1] - bg.crop[4]):(dim(timg)[1] - bg.crop[3]), bg.crop[1]:bg.crop[2]]
  timg=as.cimg(timg)
  first_Day_jpg=Calc_COG_from_cimg(timg)
  
  timg2=as.raster(load.image(file = AllFiles[length(AllFiles)]))
  
  timg2=timg2[(dim(timg2)[1] - bg.crop[4]):(dim(timg2)[1] - bg.crop[3]), bg.crop[1]:bg.crop[2]]
  timg2=as.cimg(timg2)
  Last_Day_jpg=Calc_COG_from_cimg(timg2)
  
  #load MXnet classification model
  model_loaded = mx.model.load(mx.model, iter.mx)
  
  #Cycle thru all images in the directory.
  opening_spikelets=NULL
  for(i in 1:length(AllFiles)){
    #crop image
    adjust_pic=ReadOneDayPanicle(tfile = AllFiles[i], bg.crop = bg.crop, x.move = 0, y.move = 0)
    #x
    x=round((dim(adjust_pic)[2]/44)*resize_factor, digits = 0)
    #y
    #resize image
    y=round((dim(adjust_pic)[1]/44)*resize_factor, digits = 0)
    adjust_pic=as.cimg(adjust_pic)
    adjust_pic=resize(adjust_pic, size_x = x*44, size_y = y*44)
    adjust_pic=aperm(adjust_pic, c(1,2,4,3))
    
    #structure images into matrix that will be loaded into model, by cuttung the image into 44x44 images.
    list.test=list()
    for (i in 1:x){
      for (j in 1:y){
        temp=adjust_pic[(((i-1)*44+1):(i*44)),(((j-1)*44+1):(j*44)),,]
        #plot(as.raster(as.cimg(temp)))
        list.test=c(list.test, list(temp))
      }
    }
    test.array=array(0,c(44,44,3,x*y))
    for(i in 1:(x*y)){
      test.array[,,,i]=list.test[[i]]
    }
    
    #predict open spikelets
    system.time(preds <- predict(model_loaded, test.array))
    pred.label=apply(preds,2,function(x) ifelse(x[1] > .5, yes = 0, no = ifelse(x[3] > .5, yes = 2, no=1)))
    
    #Structure the results into a table that gives open, closed, and background values for each panicle image.
    background=length(pred.label[names(table(pred.label))[1]==pred.label])
    closed=length(pred.label[names(table(pred.label))[2]==pred.label])
    open=length(pred.label[names(table(pred.label))[3]==pred.label])
    spikes=closed+open
    if(is.na(as.numeric(table(pred.label)[3]))){open=0}
    if(is.na(as.numeric(table(pred.label)[2]))){closed=0}
    if(is.na(as.numeric(table(pred.label)[2])) & is.na(as.numeric(table(pred.label)[3]))){spikes=0}
    temp=c(background,spikes,closed,open)
    print(temp)
    opening_spikelets=cbind(opening_spikelets, temp)
    
  }
  #return results
  rownames(opening_spikelets)=c("background", "spikes", "closed", "open")
  return(opening_spikelets)
  
}

#Read image and crop the image function. Day time image do not require a x.move or y.move. they will be equal to 0
ReadOneDayPanicle<-function(tfile, bg.crop, x.move, y.move)
{
  require(jpeg)
  require(imager)
  
  print(paste("Reading File: ",tfile))
  
  #read the image
  Next_img=as.raster(load.image(file = tfile))
  #crop the image
  Next_img=Next_img[((dim(Next_img)[1] - bg.crop[4])+y.move):((dim(Next_img)[1] - bg.crop[3])+y.move), (bg.crop[1]+x.move):(bg.crop[2]+x.move)]
  
  return(Next_img)
}



