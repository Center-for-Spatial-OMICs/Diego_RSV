library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)


sobj_list <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult12_Neonate12_seuobj_list.rds")


obj <- sobj_list$Adult_1
obj <- obj %>% NormalizeData()

Sobj_sub <- subset(obj, TMA_2 %in% "Ap16-3894")

cols = brewer.pal(10, 'Paired')
cols <- cols[seq_along(unique(unsup$clust))]


plt <- ImageFeaturePlot(Sobj_sub, 
                        border.color = 'white', 
                        size = 0.5,
                        features = c(top_5_auc_markers[1])) +
  #coord_flip() +
  ggtitle("") +
  theme(aspect.ratio = 1)
  

print(plt)

print(plt)




library(magick)
library(ggplot2)
library(grid)



## Create images as an example 
library(magick)

# Function to create a scatter dot image
create_scatter_image <- function(filename, shift_x = 0, shift_y = 0) {
  img <- image_blank(300, 300, "black") %>%
    image_draw()
  
  set.seed(123)  # For reproducibility
  
  # Generate random scatter points for three colors
  points(runif(30, 50, 250) + shift_x, runif(30, 50, 250) + shift_y, col = "red", pch = 16, cex = 2)
  points(runif(30, 50, 250) + shift_x, runif(30, 50, 250) + shift_y, col = "green", pch = 16, cex = 2)
  points(runif(30, 50, 250) + shift_x, runif(30, 50, 250) + shift_y, col = "blue", pch = 16, cex = 2)
  
  dev.off()
  
  image_write(img, filename)
}

# Create two images: img1 (original), img2 (shifted slightly)
create_scatter_image("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2/tmp/image1.png")
create_scatter_image("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2/tmp/image2.png", shift_x = 10, shift_y = 10)

# Display the created images
img1 <- image_read("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2/tmp/image1.png")
img2 <- image_read("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2/tmp/image2.png")

img1
img2


# Resize both images to 200x200 pixels
img1_resized <- image_scale(img1, "200x200")
img2_resized <- image_scale(img2, "200x200")

# Make img2 semi-transparent
img2_transparent <- image_colorize(img2_resized, opacity = 50, color = "white")

# Overlay the images
result <- image_composite(img1_resized, img2_transparent, operator = "atop")

# Display the result
print(result)

















setwd("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2/tmp")
library(magick)

# Function to create a scatter dot image
create_scatter_image <- function(filename, shift_x = 0, shift_y = 0) {
  img <- image_blank(300, 300, "white") %>%
    image_draw()
  
  set.seed(123)  # For reproducibility
  
  # Generate random scatter points for three colors
  points(runif(30, 50, 250) + shift_x, runif(30, 50, 250) + shift_y, col = "red", pch = 16, cex = 2)
  points(runif(30, 50, 250) + shift_x, runif(30, 50, 250) + shift_y, col = "green", pch = 16, cex = 2)
  points(runif(30, 50, 250) + shift_x, runif(30, 50, 250) + shift_y, col = "blue", pch = 16, cex = 2)
  
  dev.off()
  
  image_write(img, filename)
}

# Create aligned images
create_scatter_image("aligned_1.png")
create_scatter_image("aligned_2.png")

# Create non-aligned images (shifted by 20 pixels)
create_scatter_image("non_aligned_1.png")
create_scatter_image("non_aligned_2.png", shift_x = 20, shift_y = 20)

# Read the images
aligned_1 <- image_read("aligned_1.png")
aligned_2 <- image_read("aligned_2.png")
non_aligned_1 <- image_read("non_aligned_1.png")
non_aligned_2 <- image_read("non_aligned_2.png")

# Make the second images semi-transparent
aligned_2_transparent <- image_fx(aligned_2, expression = "a*0.5", channel = "alpha")
non_aligned_2_transparent <- image_fx(non_aligned_2, expression = "a*0.5", channel = "alpha")

# Overlay the images
aligned_overlay <- image_composite(aligned_1, aligned_2_transparent, operator = "atop")
non_aligned_overlay <- image_composite(non_aligned_1, non_aligned_2_transparent, operator = "atop")

# Display the results
print(aligned_overlay)
print(non_aligned_overlay)

# Save the results
image_write(aligned_overlay, "aligned_overlay.png")
image_write(non_aligned_overlay, "non_aligned_overlay.png")

# Create a side-by-side comparison
comparison <- image_append(c(aligned_overlay, non_aligned_overlay))
print(comparison)
image_write(comparison, "comparison.png")


setwd("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2/tmp")
ROI_img <- image_read("ROI.png")
FOV_img <- image_read("FOV.png")


ROI_img_transparent <- image_fx(ROI_img, expression = "a*0.5", channel = "alpha")


FOV_overlay <- image_composite(FOV_img, ROI_img_transparent, operator = "atop", offset = "-40-0"); print(FOV_overlay)
## +right, +down; -left,-up; units = pixels
image_write(FOV_overlay, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2/tmp/Adult1_ROI_FOV_overlay.png")


library(magick)
library(imager)
# Convert to grayscale
img <- magick2cimg(FOV_overlay)
gray_img <- grayscale(img)

# Detect white lines (grid)
edges <- cannyEdges(gray_img)
white_lines <- threshold(edges, thr = "auto")

# Find squares formed by white lines
closed <- closing(white_lines, kern = makeBrush(5, shape = "box"))
labeled <- label(closed)

# Detect greenish squares within the white-lined squares
green_channel <- channel(img, 2)
green_mask <- threshold(green_channel, thr = "auto")

# Combine labeled squares with green mask
result <- labeled * green_mask

# Visualize the results
par(mfrow = c(2, 2))
plot(img, main = "Original Image")
plot(white_lines, main = "White Lines")
plot(labeled, main = "Labeled Squares")
plot(result, main = "Final Result")




library(magick)
library(imager)
# Convert to HSV (helps with color detection)
FOV_overlay_cimg <- magick2cimg(FOV_overlay)
hsv_img <- RGBtoHSV(FOV_overlay_cimg)

# Extract channels
hue <- as.array(hsv_img[,,1,1])
saturation <- as.array(hsv_img[,,1,2])
value <- as.array(hsv_img[,,1,3])

# Identify green areas (tune thresholds if needed)
green_mask <- (hue > 0.25 & hue < 0.45) & (saturation > 0.3) & (value > 0.2)

# Identify white grid areas (tune thresholds if needed)
white_mask <- (saturation < 0.1) & (value > 0.8)

# Calculate areas
total_white_area <- sum(white_mask)
total_green_area <- sum(green_mask)

# Proportion of green inside the grid
proportion_green_in_grid <- total_green_area / total_white_area





### Registration ------
# 2. Create a binary mask of the structures in the overlay image
# This step depends on the nature of your structures. Here's a simple threshold example:
mask <- image_threshold(ROI_img, "white", "35%")

# 3. Apply the mask to the base image
result <- image_composite(FOV_img, mask, operator = "multiply")

# 4. Optionally, you can color the masked areas
colored_mask <- image_colorize(mask, opacity = 50, color = "red")
result_colored <- image_composite(FOV_img, colored_mask, operator = "over")

# 5. Display and save results
print(result)
print(result_colored)



### Removing text elements from ROI --------
# Create a mask of the text area (you may need to adjust this based on your image)
mask <- image_threshold(ROI_img, "white", "50%")

# Invert the mask
mask_inv <- image_negate(mask)

# Apply the mask to remove text
img_no_text <- image_composite(ROI_img, mask_inv, operator = "copyopacity")

# Fill in the removed area (using a simple blur method)
img_filled <- image_composite(
  img_no_text,
  image_blur(img_no_text, radius = 10, sigma = 5),
  operator = "over"
)

# Display the result
print(img_filled)




library(magick)


# Convert to grayscale
img_gray <- image_convert(ROI_img, colorspace = "gray")

# Apply edge detection to highlight text-like features
img_edges <- image_convolve(img_gray, matrix("sobel"))

# Threshold the edge image to create a mask
mask <- image_threshold(img_edges, "white", "15%")

# Invert the mask
mask_inv <- image_negate(mask)

# Dilate the mask to expand text regions
mask_dilated <- image_morphology(mask_inv, "dilate", "diamond:2")

# Apply the mask to the original image
img_no_text <- image_composite(ROI_img, mask_dilated, operator = "copyopacity")

# Fill in the removed areas
img_filled <- image_fill(img_no_text, "content-aware")

# Display and save the result
print(img_filled)
























# install.packages("RNiftyReg")
library(RNiftyReg)
# Assuming you've already performed registration
result <- niftyreg(FOV_img, ROI_img)

# Function to transform coordinates
transform_coordinates <- function(x, y, transform) {
  coords <- matrix(c(x, y, 1), ncol = 3)
  transformed <- coords %*% t(transform)
  return(list(x = transformed[1], y = transformed[2]))
}

# Example coordinates of a structure in the overlay image
overlay_coords <- list(x = 100, y = 150)

# Transform coordinates to base image space
base_coords <- transform_coordinates(overlay_coords$x, overlay_coords$y, result$forward)







ROI_img_scaled <- image_scale(ROI_img, "50%")  # Scale to 50% of original size
ROI_img_transparent <- image_fx(ROI_img_scaled, expression = "a*0.5", channel = "alpha")
FOV_overlay <- image_composite(FOV_img, ROI_img_transparent, operator = "atop", offset = "-100+0")


                      

