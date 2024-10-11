import os
from PIL import Image, ImageDraw, ImageFont

# Specify the folder containing your PNG images
folder_path = 'figs'  # Replace with your folder path

# Create a list to hold the image objects
images = []

# Loop through the expected filenames
for i in range(0, 101):  # From 0 to 100 (assuming you want 0.00 to 1.00)
    filename = f"{i/100:.2f}.png"  # Format as 0.00, 0.01, ..., 1.00
    file_path = os.path.join(folder_path, filename)
    
    # Check if the file exists
    if os.path.isfile(file_path):
        # Open the image
        img = Image.open(file_path)
        
        # Create a draw object
        draw = ImageDraw.Draw(img)
        
        # Prepare the text
        bias_text = f"Bias={filename[:-4]}"  # Remove the '.png' part
        
        # Define a font (you can specify a path to a .ttf file if needed)
        font_path = r'C:\Windows\Fonts\times.ttf'  # Adjust this to your font file path
        font_size = 32  # Change this value for a larger or smaller font
        font = ImageFont.truetype(font_path, font_size)
        
        # Get the bounding box of the text
        text_bbox = draw.textbbox((0, 0), bias_text, font)
        text_width = text_bbox[2] - text_bbox[0]
        text_height = text_bbox[3] - text_bbox[1]
        
        # Define text position (bottom-right corner)
        position = (10, img.height - text_height - 30)  # 10 pixels from the edge
        
        # Draw the text on the image
        print(bias_text)
        draw.text(position, bias_text, fill="black", font=font)

        # Append the modified image to the list
        images.append(img)

# Save the images as a GIF
if images:
    images[0].save('figs/output.gif', save_all=True, append_images=images[1:], duration=500, loop=0)
