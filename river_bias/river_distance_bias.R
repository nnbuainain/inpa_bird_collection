## Written by Nelson Buainain, 2021

## This script tests if there is a river-biased in bird sampling
#at INPA's genetic resource collection

" For that I filtered duplicate points, and points that are 40 Km close 
to each other to eliminate spatial correlation betweem them, and calculated
the distance from each point to the closest river. I then generated the same
number of random points as the original points and calculated the same distance.
At last, I generated density plots to compared the original and simulated
data and performed hypothesis tests to see if there were differences"

# Load packages, if you don't have them installed, please intall them first
#install.packages(c('sp','dismo','plyr','rgeos','ggplot2'))
library('sp')
library('dismo')
library('plyr')
library('rgeos')
library('ggplot2')


#Read table
samples = read.csv('river_bias/samples_amazonia.csv',sep=',')
head(samples)
#transform in georeferenced points
samples = SpatialPoints(samples[,c('LON_DEC','LAT_DEC')],proj4string=CRS("+proj=longlat +datum=WGS84"))

#Eliminate duplicaes within a 40km buffer
samples_filt<-remove.duplicates(samples, zero = 40.0)

#How many row is left?

samples_filt #168

#Read shapefile

amazon = shapefile('river_bias/shape_files/amazon_shape.shp')

# Set projection equals to the points
projection(amazon)<-projection(samples_filt)

#Read river shapefile

river = shapefile('river_bias/shape_files/rivers_crop.shp')
projection(river)<-projection(samples_filt)

#plot maps and points
plot(amazon)
plot(river,add=T,col='blue')
plot(samples_filt,add=T,pch=1,cex=0.5)

#Does it look good? if yes, save filtered points as csv

write.csv(samples_filt,file='river_bias/points_buffer_40km.csv')

#Generate the random points inside the amazon, number of points will be
# the same as filtered poins
set.seed(123)
random_pts = spsample(amazon,n=168,'random')
projection(random_pts)<-projection(samples_filt)

# Plot random and original points
plot(amazon)
plot(river,add=T,lwd=0.8)
plot(samples_filt,add=T,pch=16,cex=0.75,col='#0072B2')

plot(random_pts,add=T,pch=18,cex=0.75,col='#E69F00')
legend('topleft',legend = c('Original points','Random points'),fill=c('#0072B2','#E69F00'),cex=0.6,inset=0.13)


#save random points if you would like to

write.csv(random,file='river_bias/random_points.csv')

#Save plot if you would like to

png(file="river_bias/original_random_points_hq.png",width = 600,height = 400)
plot(amazon)
plot(river,add=T,lwd=0.8,)
plot(samples_filt,add=T,pch=16,cex=1.3,col='#0072B2')
plot(random_pts,add=T,pch=18,cex=1.3,col='#E69F00')
legend('topleft',legend = c('Original points','Random points'),fill=c('#0072B2','#E69F00'),cex=1.1)
dev.off()

### Distance between points and closest river

#Calculate distance to the closest river
original_dist = apply(gDistance(samples_filt, river,byid=TRUE),2,min)
random_dist = apply(gDistance(random_pts, river,byid=TRUE),2,min)

## HYPOTHESIS TEST
# Test if there is a difference among distances of original samples
# and simulated random samples to rivers

#First see if data have normal distribution
hist(original_dist)
hist(random_dist)
shapiro.test(original_dist)
shapiro.test(random_dist)

# It seems like the data is not normally distributed, also, we have a long 
# tail with some possible outliers, wilcoxon test might be a better alternative
# to a student test

wilcox.test(original_dist, random_dist, paired = FALSE)

# The test had a very high value for W, thus different from 0 and very low
# p-value indicating that there was a difference among the two datasets

# Ploting graphs to see distribution of these distances
#create dataframes with classification, this will make easier to make
#plots with ggplot2

river_dist_original = data.frame(dist=original_dist,Points='original')
river_dist_random = data.frame(dist=random_dist,Points='random')

#bind the two df into one

river_dist = rbind(river_dist_original,river_dist_random)
river_dist

ggplot(river_dist, aes(x = dist, fill = Points)) + geom_density(alpha = 0.5)

# Calculate median of the two datasets
median <- ddply(river_dist, "Points", summarise, dist.median=median(dist))
head(median)

ggplot(river_dist, aes(x = dist, fill = Points)) +
  geom_density(alpha = 0.5) +
  geom_vline(data=median, aes(xintercept=dist.median, color=Points),
             linetype="dashed")+
  scale_color_manual(values=c("#0072B2", "#E69F00"))+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  labs(title="River biased sampling",x="Distance to the closest river", y = "Density")+
annotate("text",label="Wilcoxon test – p<0.001",x=3,y=1,lwd=6)

# Ok this is good but the shapefile we used was too simplistic because
# it only has the main river courses, let's try using a more detailed 
# and refined river shapefile

###### Try a different river shape

river_fine = shapefile('river_bias/shape_files/rivers_refined.shp')
projection(river_fine)<-projection(samples_filt)

plot(amazon)
plot(river_fine,add=T)

### Distance between points and closest river

#Calculate distance to the closes river

original_dist_fine = apply(gDistance(samples_filt, river_fine,byid=TRUE),2,min)
random_dist_fine = apply(gDistance(random_pts, river_fine,byid=TRUE),2,min)

## HYPOTHESIS TEST

hist(original_dist_fine)
shapiro.test(original_dist_fine)
hist(random_dist_fine)
shapiro.test(random_dist_fine)

wilcox.test(original_dist_fine, random_dist_fine, paired = FALSE)

#create dataframes with classification
river_dist_original_fine = data.frame(dist=original_dist_fine,Points='original')
river_dist_random_fine = data.frame(dist=random_dist_fine,Points='random')

#bind the two df

river_dist_fine = rbind(river_dist_original_fine,river_dist_random_fine)
river_dist_fine


median_fine <- ddply(river_dist_fine, "Points", summarise, dist.median=median(dist))
head(median_fine)

ggplot(river_dist_fine, aes(x = dist, fill = Points)) +
  geom_density(alpha = 0.5) +
  geom_vline(data=median_fine, aes(xintercept=dist.median, color=Points),
             linetype="dashed")+
  scale_color_manual(values=c("#0072B2", "#E69F00"))+
  scale_fill_manual(values=c("#0072B2", "#E69F00"))+
  labs(title="River biased sampling",x="Distance to the closest river (in º)", y = "Density")+
  annotate("text",label="Wilcoxon test – p<0.001",x=0.8,y=4.5,lwd=6)+
  theme(legend.position = c(0.8,0.8))

png(file="river_bias/density_plot_river_distances.png",width = 400,height = 400)
dev.off()
