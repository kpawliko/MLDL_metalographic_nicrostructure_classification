import os
import random
import cv2
import numpy as np



#running RCALib testing tool
def runRCA(variant, neighbourhood, grainsCount, outputName):
    os.system(f'../build/benchmark \
        floating-percision=1 \
        material-version=state-changed-flag \
        cell-generator=constant \
            x-min=-1.0 \
            x-max=1.0 \
            y-min=-1.0 \
            y-max=1.0 \
            cell-count=200000 \
        grain-arrangement=random \
        variant={variant}\
        grain-count={grainsCount} \
        neighbourhood={neighbourhood} \
            radius=0.03 \
            ra=0.03\
            rb=0.015\
        rule=grain-growth \
        output=images \
            image-width=300 \
            image-height=300 \
            images-directory={outputName} \
        method=FixedGrid \
            space-redundancy=0.2')



def detectBoundaries(input_path, outputPath):
    
    img = cv2.imread(input_path) 
    img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    kernel = np.ones((5, 5), np.uint8)
    img_dilation = cv2.dilate(img_gray, kernel, iterations=1) 
    img_blur = cv2.GaussianBlur(img_dilation, (3,3), 0) 
    edges = cv2.Canny(image=img_blur, threshold1=1, threshold2=20) 	
    cv2.imwrite(outputPath, edges) 




variantsList = ["most", "nearest", "first"]
neighbourhoodsList = ["circle", "ellipse"]

# runRCA("first", "circle", 400, "t1")
# runRCA("nearest", "circle", 400, "t2")

# for i in range(50):
#     runRCA("first", "circle", 170, f'./structs/a{i}')
#     detectBoundaries(f'./structs/a{i}/final.bmp', f'train_data/desired_structs/a{i}.bmp')





for i in range(20):
    runRCA(random.choice(variantsList), random.choice(neighbourhoodsList), 2000, f'./structs2/b{i}')
    detectBoundaries(f'./structs2/b{i}/final.bmp', f'train_data/other_structs/b{i}.bmp')

# for i in range(20):
#     runRCA(random.choice(variantsList), random.choice(neighbourhoodsList), random.randint(350, 450), f'./structs2/a{i+20}')
#     detectBoundaries(f'./structs2/a{i+20}/final.bmp', f'train_data/other_structs/a{i+20}.bmp')

# for i in range(20):
#     runRCA(random.choice(variantsList), random.choice(neighbourhoodsList), i+10, f'./structs2/a{i+20}')
#     detectBoundaries(f'./structs2/a{i+20}/final.bmp', f'train_data/other_structs/a{i+20}.bmp')





# Importing all necessary libraries
