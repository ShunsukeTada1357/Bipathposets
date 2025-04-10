#In this file, we give a run test for image data. 
#We explain how  
#(1) to load image data from a folder.
#(2) to change image data into 0-1 matrices.
#Finally, we compute bipath persistence diagrams of bipath filtrations of image data.
##################################################
#(1) how to read image files automatically from a folder.
using Glob 
# Read all PNG files from two folders.
# For example...
img_dir_up = raw"C:\Users\shunsuke\.julia\packages\Bipathposets\zo6pb\test\imagesUp15/"
img_dir_down = raw"C:\Users\shunsuke\.julia\packages\Bipathposets\zo6pb\test\imagesDown15/"
#↑ change to appropriate paths.

img_paths_up = glob("*.PNG", img_dir_up)
img_paths_down = glob("*.PNG", img_dir_down)
####################################################
#(2)
function image_to_binary_matrix(mat; threshold=0.5) #mat is a matrix (a_ij) 0≤ a_ij ≤ 1. it represents gray image.   
    # 各ピクセルの値が threshold より大きければ 1, そうでなければ 0 を返す
    return map(x -> x.val > threshold ? 1 : 0, mat)
end

"""
Example
mat = Gray.(load("1_15.png")) # image data -> matrix gray scale
WhiteBlackmat = image_to_binary_matrix(mat)  # 0-1 matrix
display(Gray.(mat))
"""

#Finally, we compute bipath PD. 
using Images
using ImageView
# 画像を読み込こむ → gray scale -> 0-1 matrix (配列に格納)
#imgs = [load(path) for path in img_paths]
imgsUp = [image_to_binary_matrix(Gray.(load(path))) for path in img_paths_up]
imgsDown = [image_to_binary_matrix(Gray.(load(path))) for path in img_paths_down]

# Compute bipath persistence diagram!
aa = Bipath_image(imgsUp,imgsDown)
Bipathposets.plotintlist(aa,0)
Bipathposets.plotintlist(aa,1)