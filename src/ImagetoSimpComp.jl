using Combinatorics
#using Images
#using ImageView
########################################################
## 0-1 matrix → simplicial complex
## (Future work: Image data = matrix → cubical complex for fast computation) 
## We consider image data as 0-1 matrix.

function ImagetoFSC(mats...)  #Input are matrices, whose size are the same.
    row, col = size(mats[1])  # row: 行数、 col: 列数
    S = []  # output list
    VE = [] # list for processd.
    n = length(mats)
    for l in 1:n
       for i in 1:row
           for j in 1:col
               if mats[l][i, j] == 0  # 黒ピクセルの場合
                   # 4頂点の番号を計算
                   a_val = (col+1)*(i-1) + j          #a-b
                                                      #|\|
                   b_val = (col+1)*(i-1) + j + 1      #c-d
                   c_val = (col+1)*i + j
                   d_val = (col+1)*i + j + 1

                   #Get simplicial complex consisting of {{a},{b},{c},{d},{a,b},{a,c},{a,d},{b,d},{c,d}}
                   #(Future work: it should be cubical complex given by {a,b,c,d}.
                   P = collect(powerset([a_val, b_val, c_val, d_val],1,3))
                   setdiff!(P, [[b_val,c_val], [a_val, b_val,c_val],[ b_val,c_val, d_val] ] )
                   
                   for subset in P
                       if !(subset in VE)
                           append!(S,[[vcat(subset...),l]])
                           push!(VE, subset)
                        
                       end
                   end
               end
           end
       end
    end
    return S
end
"""
#Example
#Let A and B be the image data given by 
A = [0 1 1; 0 1 0; 0 0 0]
B = [0 0 0; 0 1 0; 0 0 0]
#Then, we obtain
aa=ImagetoFSC(A,B)
println(aa)
"""

#Take union and intersection of images.
function combine_images(Img...) #Inputs are matrices of the same size
    r, c = size(Img[1])
    Union = zeros(Int, r, c)
    Intersection = zeros(Int, r, c)
    n = length(Img)
    for k in 1:n-1
    for i in 1:r
        for j in 1:c
            if Img[k][i, j] == 0 && Img[k+1][i, j] == 0
                Union[i, j] = 0
                Intersection[i, j] = 0
            elseif Img[k][i, j] == 0 || Img[k+1][i, j] == 0
                Intersection[i, j] = 1
                Union[i, j] = 0
            else
                Union[i, j] = 1
                Intersection[i, j] = 1
            end
        end
    end
    end
    return Union, Intersection
end
"""
#Example usinf 20*20 pixel images.
mat = zeros(Int, 20, 20)
mat[5:10, 5:10] .= 1  # 中央の領域を1（白）に設定
mat2 = zeros(Int, 20, 20)
mat2[7:12, 7:12] .= 1  # 中央の領域を1（白）に設定
mat[1][1]
aa = combine_images(mat,mat2,mat)
display(Gray.(aa[1]))
display(Gray.(aa[2]))
"""

#Make bipath filtration of image data.
function SemiFImgtoFImg(Fa,Fb) #Fa matrices, Fb matrices
    X, Y = copy(Fa), copy(Fb)
    X[1] = combine_images(X[1], Y[1])[2] # take intersection 
    Y[1] = X[1]
    X[length(X)] = combine_images(last(X), last(Y))[1] # take unioun
    Y[length(Y)] = X[length(X)] 
    return X, Y
end
"""
Example
mat = zeros(Int, 20, 20)
mat[5:10, 5:10] .= 1  # 中央の領域を1（白）に設定
mat2 = zeros(Int, 20, 20)
mat2[7:12, 7:12] .= 1  # 中央の領域を1（白）に設定
mat[1][1]
aa = combine_images(mat,mat2,mat)

Fa = [aa[1], mat, aa[2] ]
Fb = [aa[1], mat2, aa[2] ]
"""

function Bipath_image(MatsUp,MatsDown)
    X = SemiFImgtoFImg(MatsUp,MatsDown)[1]
    Y =SemiFImgtoFImg(MatsUp,MatsDown)[2]
    FSCa = [ImagetoFSC(X...) , length(MatsUp)]
    FSCb = [ImagetoFSC(Y...) , length(MatsDown) ]
    return Bipathposets.interval_decomposition(FSCa,FSCb)
end


"""
#Example1
mat = zeros(Int, 20, 20)
mat[5:10, 5:10] .= 1  # 中央の領域を1（白）に設定
mat2 = zeros(Int, 20, 20)
mat2[7:12, 7:12] .= 1  # 中央の領域を1（白）に設定
mat[1][1]
aa = combine_images(mat,mat2,mat)
Fa = [aa[1], mat, aa[2] ]
Fb = [aa[1], mat2, aa[2] ]
bb = Bipath_image(Fa,Fb)
Bipathposets.plotintlist(bb,0)
Bipathposets.plotintlist(bb,1)
"""



###Below is not needed now.
"""
function FImgstoFimg(Fimgs...)
    Fimg = []
    n = length(Fimgs)
    for i in 1:n
        A = [Fimgs[k][i] for k in 1:n]
        append!(Fimg, A)
    end
    return combine_images(Fimg...)
end


function BFImgstoBFImg(BipathFImgs)
    n = length(BipathFImgs) 
    Fimgsa = [BipathFImgs[i][1] for i in 1:n]
    Fimgsb = [BipathFImgs[i][2] for i in 1:n]

    Fimgsa = FImgstoFimg(Fimgsa...)
    Fimgsb = FImgstoFimg(Fimgsb...)

    return SemiFImgtoFImg(Fimgsa,Fimgsb)
end

"""