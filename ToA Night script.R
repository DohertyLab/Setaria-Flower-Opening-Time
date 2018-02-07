##How to calculate daytime image flower opening time##
# InDir is a directory of all the panicle images for one day
# By is the sequence of images to take from the directory with images (ie By=30 would take every 30th image in the directory)
# mx.model is the MXnet model used for classification
# last_light_jpg T or F if the directory includes a picture of the daytime image to isolate spikelet locations
# iter.mx is the iteriation of the MXnet model used
# resize_factor is the numeric between 0-1 that reduces the size of the image by the factor.

timelapse_night_spikes=ReadMulti_NightPanicle(InDir = TheDir,By = 1, 
                                               last_light_jpg = T,
                                               mx.model = '/mnt/scratch/jigar/flower_opening_learning/Inception_Small-open_closed', iter.mx = 85, resize_factor = .7)

ReadMulti_NightPanicle<-function(InDir, By=1, last_light_jpg=T, mx.model=NULL, iter.mx=NULL,  resize_factor=.7)
{
  
  #read images in directory
  AllFiles=dir(path = InDir, pattern = ".jpg$", full.names = T)
  AllFiles<-sort(AllFiles)
  AllFiles<-AllFiles[seq(1,length(AllFiles),by = By)]
  names(AllFiles)<-AllFiles
  
  timg=as.raster(load.image(file = AllFiles[2]))
  plot(timg)
  
  #crop images and adjust by change in center of mass between images
  #calc_COG_from_cimg function calcuations the change in center of mass between images
  #ReadOneNightPanicle2 reads in an image and adjusts the image by center of mass and crops the image by the selected area
  bg.crop = base::as.vector(extent(select(raster(x = AllFiles[2], band = 2))))
  timg=timg[(dim(timg)[1] - bg.crop[4]):(dim(timg)[1] - bg.crop[3]), bg.crop[1]:bg.crop[2]]
  timg=as.cimg(timg)
  first_night_jpg=Calc_COG_from_cimg(timg)
  
  timg2=as.raster(load.image(file = AllFiles[length(AllFiles)]))
  
  timg2=timg2[(dim(timg2)[1] - bg.crop[4]):(dim(timg2)[1] - bg.crop[3]), bg.crop[1]:bg.crop[2]]
  timg2=as.cimg(timg2)
  Last_night_jpg=Calc_COG_from_cimg(timg2)
  
  x.move=(Last_night_jpg[1]-first_night_jpg[1])
  y.move=(Last_night_jpg[3]-first_night_jpg[3])
  
  x.move=seq(0, x.move, length.out = length(AllFiles))
  y.move=seq(0, y.move, length.out = length(AllFiles))
  
  
  for(i in 1:length(AllFiles)){
    adjust_pic=ReadOneNightPanicle2(tfile = AllFiles[i], bg.crop = bg.crop, x.move = x.move[i], y.move = y.move[i])
    save(adjust_pic, file =  paste(names(AllFiles)[i],"adjust.rda", sep = "_"))
  }
  
  # Uses the last day image to find the location of spikelets.
  if((last_light_jpg)==T){
    
    model_loaded = mx.model.load(mx.model, iter.mx)
    
    last_light_img=as.raster(load.image(file = AllFiles[1]))
    last_light_img=last_light_img[(dim(last_light_img)[1] - bg.crop[4]):(dim(last_light_img)[1] - bg.crop[3]), bg.crop[1]:bg.crop[2]]
    last_light_img=as.cimg(last_light_img)
    #x
    x=round((dim(last_light_img)[2]/44)*resize_factor, digits = 0)
    #y
    y=round((dim(last_light_img)[1]/44)*resize_factor, digits = 0)
    
    last_light_img=resize(last_light_img, size_x = x*44, size_y = y*44)
    last_light_img=aperm(last_light_img, c(1,2,4,3))
    last_light_img.list=list()
    for (i in 1:x){
      for (j in 1:y){
        temp=last_light_img[(((i-1)*44+1):(i*44)),(((j-1)*44+1):(j*44)),,]
        #plot(as.raster(as.cimg(temp)))
        last_light_img.list=c(last_light_img.list, list(temp))
        rm(temp)
      }
    }
    img.array=array(0,c(44,44,3,x*y))
    print(x)
    print(y)
    for(i in 1:(x*y)){
      img.array[,,,i]=last_light_img.list[[i]]
    }
    
    #saves the location of where the spikelets are in the image
    preds <- predict(model_loaded, img.array)
    pred.label=apply(preds,2,function(x) ifelse(x[1] > .5, yes = 0, no = ifelse(x[3] > .5, yes = 2, no=1)))
    print(table(pred.label))
    
    
    # loop crops out each spikelet from the time series. And puts each spikelet into a its own series of images. 
    timelapse_night_spikes=list()
    for(a in 1:length(list.files(InDir, pattern = "_adjust.rda"))){
      
      print(paste("timepoint subset #", a, sep = ""))
      
      night_img=load(list.files(InDir, pattern = "_adjust.rda", full.names = T)[a])
      night_img=get(night_img)
      night_img=as.cimg((night_img))
      night_img<-grayscale(night_img)[,,1,1]
      #x
      x=round((dim(night_img)[1]/44)*resize_factor, digits = 0)
      #y
      y=round((dim(night_img)[2]/44)*resize_factor, digits = 0)
      
      night_img=resize(as.cimg(night_img), size_x = x*44, size_y = y*44)  
      
      #save each spikelet timeseries in a list. each object in the list is a time series of one spikelet
      night_img.list=list()
      for (i in 1:x){
        for (j in 1:y){
          temp=night_img[(((i-1)*44+1):(i*44)),(((j-1)*44+1):(j*44)),,]
          night_img.list=c(night_img.list, list(temp))
          rm(temp)
        }
      }
      night_img.list=night_img.list[which(pred.label==1)]
      night_img.list=lapply(night_img.list, function(x) as.raster(as.cimg(x)))
      timelapse_night_spikes[[a]]<-night_img.list
    }
    
    # names and saves the spikelet timeseries list
    names(timelapse_night_spikes)<-list.files(InDir, pattern = "_adjust.rda")
    i <- 1:length(timelapse_night_spikes)
    j <- 1:length(timelapse_night_spikes[[1]])
    swap<-lapply(j, function(j) lapply(i, function(i) timelapse_night_spikes[[i]][[j]]))
    swap=lapply(swap,function(x) setNames(x, names(timelapse_night_spikes)))
    timelapse_night_spikes=swap
    rm(swap)
    save(timelapse_night_spikes, file = "timelapse_night_spikes.rda")
    return(timelapse_night_spikes)
    save(pred.label, file="pred.label.rda")
  }
  
}

# helper fucntion for nightime
#calc_COG_from_cimg function calcuations the change in center of mass between images
Calc_COG_from_cimg<-function(cimg){
  cimg<-grayscale(cimg)[,,1,1]
  cimg=as.asc(cimg)
  g=COGravity(cimg)
  return(g)
}

#ReadOneNightPanicle2 reads in an image and adjusts the image by center of mass and crops the image by the selected area
ReadOneNightPanicle2<-function(tfile, bg.crop, x.move, y.move)
{
  require(jpeg)
  require(imager)
  
  print(paste("Reading File: ",tfile))
  
  print(paste("x move", x.move))
  print(paste("y move", y.move))
  
  Next_img=as.raster(load.image(file = tfile))
  Next_img=Next_img[((dim(Next_img)[1] - bg.crop[4])+y.move):((dim(Next_img)[1] - bg.crop[3])+y.move), (bg.crop[1]+x.move):(bg.crop[2]+x.move)]
  
  return(Next_img)
}


##uses the output of ReadMulti_NightPanicle to calculate the find spikelets that open during the night. 
#timelapse_night_spikes_list is list where each each object in the list is timeseries of images of a spikelet at night
#pics_length_to_mean is the length of spikelets to average. (ie if equal to 5 every 5 images will be averaged.)
# base_to_test_interval is the number of images that will be compared to an initial timepoint. for example this is equal to 30, then each images will be bloack off into groups of 30 and each image in that block will be substracted from the first image in that block
Night_basleine_average2<-function(timelapse_night_spikes_list, pics_length_to_mean, base_to_test_interval){
  require(doParallel)  
  require(foreach)
  percent_change=NULL
  for(i in seq(1,length(timelapse_night_spikes_list),base_to_test_interval)){
    print(paste("next set", i))
    img1=calc_mean_of_time.spikes(list_of_spikes_to_mean =timelapse_night_spikes_list[i:(i+pics_length_to_mean)])
    img1=as.cimg(img1)
    img1=resize(img1,size_x = 480, size_y = 480)
    
    test=NULL
    for( j in seq(1,base_to_test_interval,pics_length_to_mean)){

      
      if((i+j+pics_length_to_mean)<length(timelapse_night_spikes_list)){
        img2=calc_mean_of_time.spikes(list_of_spikes_to_mean = timelapse_night_spikes_list[(i+(j-1)):((i+(j-1))+pics_length_to_mean)])
        img2=as.cimg(img2)
        img2=resize(img2,size_x = 480, size_y = 480)
        temp=(sum(img2+128 <= img1) + sum(img2-128 >= img1))/length(img1)
      }
      test=c(test,temp)
    }
    print(test)
    percent_change=c(percent_change,test)
  }
  return(percent_change)
}

##make GIF of timelapse_night_spikes
#uses the output of ReadMulti_NightPanicle to make a gif of each individual spikelets
system(paste("mkdir", paste(TheDir,"/gif", sep = ""), sep = " "))
setwd(paste(TheDir,"/gif", sep = ""))

i <- 1:length(timelapse_night_spikes)
j <- 1:length(timelapse_night_spikes[[1]])
swap<-lapply(j, function(j) lapply(i, function(i) timelapse_night_spikes[[i]][[j]]))
swap=lapply(swap,function(x) setNames(x, names(timelapse_night_spikes)))

cl <- makeCluster(10)  
registerDoParallel(cl)
foreach(i=c(1:length(swap[[1]])), .packages=c('imager','raster','XML','Cairo')) %dopar% {
  dir.create(paste("spike",i, sep = ""))
  setwd(paste("spike",i, sep = ""))
  for(j in 1:length(swap)){
    CairoPNG(filename = paste(list.files(path = TheDir, pattern = "_adjust.rda")[j],paste(i,"spike.png", sep = ""), sep = "_"))
    plot(swap[[j]][[i]])
    dev.off()
  }
  system("/usr/local/bin/mogrify -font Liberation-Sans -fill white -undercolor '#00000080' \ -pointsize 26 -gravity NorthEast -annotate +10+10 %t *.png")
  system("/usr/local/bin/convert -delay 2 *.png spike.gif")
  setwd("..")
}
stopCluster(cl)
