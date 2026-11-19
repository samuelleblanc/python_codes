#!/usr/bin/env python
# coding: utf-8

# In[2]:


import collections
import collections.abc
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.shapes import MSO_SHAPE
from pptx.dml.color import RGBColor

# Initialize slide canvas (16:9 widescreen layout)
prs = Presentation()
prs.slide_width = Inches(13.33)
prs.slide_height = Inches(7.5)
slide_layout = prs.slide_layouts[6] # Blank canvas layout
slide = prs.slides.add_slide(slide_layout)

# Define Scientific Color Scheme
SKY_BLUE = RGBColor(215, 235, 250)
PBL_GRAY = RGBColor(220, 225, 230)
TEXT_DARK = RGBColor(30, 30, 30)
ARROW_BLUE = RGBColor(70, 130, 180)
LINE_AEROSOL = RGBColor(180, 40, 40)
LINE_NO2 = RGBColor(40, 80, 180)
LINE_O3 = RGBColor(40, 140, 60)

# --- 1. Background Sky and PBL Gradients ---
sky = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.5), Inches(0.5), Inches(8.0), Inches(4.5))
sky.fill.solid()
sky.fill.fore_color.rgb = SKY_BLUE
sky.line.fill.background()

pbl = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.5), Inches(5.0), Inches(8.0), Inches(1.5))
pbl.fill.solid()
pbl.fill.fore_color.rgb = PBL_GRAY
pbl.line.fill.background()

# --- 2. Title & Key Platform Labels ---
def add_clean_label(text, left, top, width, height, size=11, bold=True):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(size)
    p.font.bold = bold
    p.font.color.rgb = TEXT_DARK
    return txBox

add_clean_label("TEMPO Geostationary Satellite", Inches(3.0), Inches(0.6), Inches(4.0), Inches(0.4), size=12)
add_clean_label("NASA P-3", Inches(1.0), Inches(3.3), Inches(1.5), Inches(0.3))
add_clean_label("NASA 777", Inches(4.5), Inches(2.2), Inches(1.5), Inches(0.3))

# --- 3. Aircraft & Satellite Interactive Vectors ---
# NASA P-3 Vertical Profiling Indication
arrow_p3 = slide.shapes.add_shape(MSO_SHAPE.UP_DOWN_ARROW, Inches(2.8), Inches(2.8), Inches(0.4), Inches(1.8))
arrow_p3.fill.solid()
arrow_p3.fill.fore_color.rgb = ARROW_BLUE
add_clean_label("P-3 Profiling Range", Inches(3.3), Inches(3.5), Inches(1.5), Inches(0.4), size=10, bold=False)

# TEMPO Downlooking Remote Sensing Vectors
arrow_sat1 = slide.shapes.add_shape(MSO_SHAPE.DOWN_ARROW, Inches(4.0), Inches(1.2), Inches(0.2), Inches(0.8))
arrow_sat2 = slide.shapes.add_shape(MSO_SHAPE.DOWN_ARROW, Inches(4.3), Inches(1.2), Inches(0.2), Inches(0.8))
for arr in [arrow_sat1, arrow_sat2]:
    arr.fill.solid()
    arr.fill.fore_color.rgb = ARROW_BLUE
    arr.line.fill.background()

# --- 4. Callout Bubbles & Measurement Panels ---
# Airborne Sunphotometer Indicator Bubble
bubble_sun = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.6), Inches(1.5), Inches(1.8), Inches(1.0))
bubble_sun.fill.solid()
bubble_sun.fill.fore_color.rgb = RGBColor(255, 255, 255)
add_clean_label("Airborne\nSunphotometer", Inches(0.7), Inches(1.6), Inches(1.6), Inches(0.8), size=10)

# Total Column Measurement Label
panel_meas = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(4.8), Inches(1.1), Inches(2.2), Inches(0.9))
panel_meas.fill.solid()
panel_meas.fill.fore_color.rgb = RGBColor(255, 255, 255)
add_clean_label("Total Column\nMeasurements:\nNO2, O3, AOD", Inches(4.9), Inches(1.15), Inches(2.0), Inches(0.8), size=10)

# --- 5. Data Assimilation Block (WRF-Chem CTM Loop) ---
# Outer Assimilation Loop Arrow
loop_arrow = slide.shapes.add_shape(MSO_SHAPE.BLOCK_ARC, Inches(9.0), Inches(1.5), Inches(3.5), Inches(3.5))
loop_arrow.fill.solid()
loop_arrow.fill.fore_color.rgb = ARROW_BLUE
loop_arrow.line.fill.background()

# Central Core Label Box
add_clean_label("WRF-Chem\nCTM", Inches(9.7), Inches(2.8), Inches(2.2), Inches(0.8), size=16, bold=True)
add_clean_label("Model-Measurement\nSynthesis", Inches(9.6), Inches(3.8), Inches(2.4), Inches(0.5), size=10, bold=False)

# --- 6. The Vertical Profiling Cartesian Axis Elements ---
def draw_profile_axis(x_start, name, line_color):
    # Base Box Frame
    frame = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x_start), Inches(1.0), Inches(1.2), Inches(5.5))
    frame.fill.solid()
    frame.fill.fore_color.rgb = RGBColor(255, 255, 255)
    frame.line.color.rgb = TEXT_DARK
    
    # Title Header
    add_clean_label(name, Inches(x_start), Inches(0.6), Inches(1.2), Inches(0.4), size=11, bold=True)
    
    # Altitude Vertical Grid lines
    add_clean_label("10", Inches(x_start - 0.3), Inches(2.0), Inches(0.3), Inches(0.3), size=9, bold=False)
    add_clean_label("5", Inches(x_start - 0.3), Inches(4.0), Inches(0.3), Inches(0.3), size=9, bold=False)
    add_clean_label("2", Inches(x_start - 0.3), Inches(5.2), Inches(0.3), Inches(0.3), size=9, bold=False)
    add_clean_label("0 (km)", Inches(x_start - 0.4), Inches(6.3), Inches(0.4), Inches(0.3), size=9, bold=False)

# Render Profiles Side-by-Side on Right Dashboard Area
draw_profile_axis(8.8, "Aerosol", LINE_AEROSOL)
draw_profile_axis(10.3, "NO2", LINE_NO2)
draw_profile_axis(11.8, "O3", LINE_O3)

# Save presentation ready to be modified directly inside PowerPoint
prs.save("Atmospheric_Remote_Sensing_Diagram.pptx")
print("Vector PowerPoint diagram generated successfully as 'Atmospheric_Remote_Sensing_Diagram.pptx'.")


# In[6]:


import os
os.curdir


# In[8]:


os.path.abspath(os.curdir)


# In[9]:


import collections
import collections.abc
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.shapes import MSO_SHAPE
from pptx.dml.color import RGBColor

# Initialize slide canvas (16:9 widescreen layout)
prs = Presentation()
prs.slide_width = Inches(13.33)
prs.slide_height = Inches(7.5)
slide_layout = prs.slide_layouts[6] # Blank slide layout
slide = prs.slides.add_slide(slide_layout)

# Define Extended Scientific Color Scheme
SKY_BLUE = RGBColor(195, 225, 245)
PBL_GRAY = RGBColor(215, 220, 225)
GROUND_DARK = RGBColor(100, 105, 110)
TEXT_DARK = RGBColor(30, 30, 30)
SUN_YELLOW = RGBColor(255, 215, 0)
ARROW_BLUE = RGBColor(65, 125, 185)
WHITE = RGBColor(255, 255, 255)
MOUNTAIN_BROWN = RGBColor(145, 130, 115)

LINE_AEROSOL = RGBColor(180, 40, 40)
LINE_NO2 = RGBColor(40, 80, 180)
LINE_O3 = RGBColor(40, 140, 60)

# --- 1. Background Layers (Sky, Planetary Boundary Layer, Ground) ---
sky = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.5), Inches(0.5), Inches(8.0), Inches(4.5))
sky.fill.solid(); sky.fill.fore_color.rgb = SKY_BLUE; sky.line.fill.background()

pbl = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.5), Inches(5.0), Inches(8.0), Inches(1.5))
pbl.fill.solid(); pbl.fill.fore_color.rgb = PBL_GRAY; pbl.line.fill.background()

ground = slide.shapes.add_shape(MSO_SHAPE.PARALLELOGRAM, Inches(0.5), Inches(6.5), Inches(8.0), Inches(0.6))
ground.fill.solid(); ground.fill.fore_color.rgb = GROUND_DARK; ground.line.fill.background()

# --- 2. Natural Features (Sun, Clouds, Background Mountains) ---
sun = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.2), Inches(0.2), Inches(1.2), Inches(1.2))
sun.fill.solid(); sun.fill.fore_color.rgb = SUN_YELLOW; sun.line.fill.background()

# Background Mountains
mtn1 = slide.shapes.add_shape(MSO_SHAPE.ISOSCELES_TRIANGLE, Inches(3.2), Inches(4.2), Inches(1.8), Inches(1.5))
mtn1.fill.solid(); mtn1.fill.fore_color.rgb = MOUNTAIN_BROWN; mtn1.line.fill.background()
mtn2 = slide.shapes.add_shape(MSO_SHAPE.ISOSCELES_TRIANGLE, Inches(4.5), Inches(4.0), Inches(2.2), Inches(1.7))
mtn2.fill.solid(); mtn2.fill.fore_color.rgb = MOUNTAIN_BROWN; mtn2.line.fill.background()

# Clouds
cloud1 = slide.shapes.add_shape(MSO_SHAPE.CLOUD, Inches(2.5), Inches(0.8), Inches(1.5), Inches(0.7))
cloud1.fill.solid(); cloud1.fill.fore_color.rgb = WHITE; cloud1.line.fill.background()

# --- 3. Urban Area Features (Buildings & Roadways) ---
def add_building(left, top, width, height, face_color):
    bldg = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(left), Inches(top), Inches(width), Inches(height))
    bldg.fill.solid(); bldg.fill.fore_color.rgb = face_color
    bldg.line.color.rgb = RGBColor(50, 50, 50)
    return bldg

# Rendering various cityscape heights
add_building(1.0, 5.5, 0.4, 1.0, RGBColor(180, 190, 200))
add_building(1.5, 5.2, 0.5, 1.3, RGBColor(160, 170, 180))
add_building(2.2, 4.8, 0.6, 1.7, RGBColor(140, 150, 160))
add_building(2.9, 4.5, 0.5, 2.0, RGBColor(190, 195, 200))

# Factory Smoke Stacks with simplified plume lines
factory = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(6.5), Inches(5.3), Inches(0.8), Inches(1.2))
factory.fill.solid(); factory.fill.fore_color.rgb = RGBColor(130, 130, 130)
plume = slide.shapes.add_shape(MSO_SHAPE.UP_ARROW, Inches(6.7), Inches(4.6), Inches(0.4), Inches(0.7))
plume.fill.solid(); plume.fill.fore_color.rgb = WHITE

# Airflow/Transport Arrow
trans_arrow = slide.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, Inches(5.8), Inches(5.1), Inches(1.5), Inches(0.2))
trans_arrow.fill.solid(); trans_arrow.fill.fore_color.rgb = TEXT_DARK

# --- 4. Platforms & Remote Sensing Instrument Assets ---
# NASA P-3 Body Shape Representative
p3_body = slide.shapes.add_shape(MSO_SHAPE.CHEVRON, Inches(1.2), Inches(2.8), Inches(1.6), Inches(0.4))
p3_body.fill.solid(); p3_body.fill.fore_color.rgb = WHITE; p3_body.line.color.rgb = ARROW_BLUE

# Direct Solar Tracking Vector
solar_beam = slide.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, Inches(1.0), Inches(1.2), Inches(1.5), Inches(1.4))
solar_beam.rotation = 45
solar_beam.fill.solid(); solar_beam.fill.fore_color.rgb = LINE_AEROSOL

# TEMPO Satellite Base Core
tempo_core = slide.shapes.add_shape(MSO_SHAPE.HEXAGON, Inches(5.2), Inches(0.8), Inches(0.8), Inches(0.8))
tempo_core.fill.solid(); tempo_core.fill.fore_color.rgb = RGBColor(90, 100, 110)
# Solar Arrays Left/Right
panel_l = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(4.2), Inches(1.0), Inches(1.0), Inches(0.4))
panel_r = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(6.0), Inches(1.0), Inches(1.0), Inches(0.4))
for panel in [panel_l, panel_r]:
    panel.fill.solid(); panel.fill.fore_color.rgb = ARROW_BLUE

# --- 5. Interactive Labels & Functional Loops ---
def add_clean_label(text, left, top, width, height, size=11, bold=True, center=False):
    txBox = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(size)
    p.font.bold = bold
    p.font.color.rgb = TEXT_DARK
    if center:
        p.alignment = 1
    return txBox

add_clean_label("TEMPO Geostationary Satellite", 4.5, 0.4, 3.5, 0.4, size=12)
add_clean_label("NASA P-3", 1.2, 3.2, 1.5, 0.3)
add_clean_label("Airflow/Transport", 5.8, 4.8, 1.5, 0.3, size=10, bold=False)
add_clean_label("Urban Area", 3.5, 6.6, 2.0, 0.4, size=12, bold=True)

# P-3 Profiling Range Dual-Arrow Paths
arrow_p3_1 = slide.shapes.add_shape(MSO_SHAPE.UP_DOWN_ARROW, Inches(3.6), Inches(2.5), Inches(0.3), Inches(1.8))
arrow_p3_2 = slide.shapes.add_shape(MSO_SHAPE.UP_DOWN_ARROW, Inches(4.1), Inches(2.7), Inches(0.2), Inches(1.4))
arrow_p3_1.fill.solid(); arrow_p3_1.fill.fore_color.rgb = ARROW_BLUE
arrow_p3_2.fill.solid(); arrow_p3_2.fill.fore_color.rgb = RGBColor(120, 160, 200)

add_clean_label("P-3 Sounding\nRange", 3.5, 1.8, 1.2, 0.6, size=10, center=True)
add_clean_label("P-3\nProfiling\nRange", 4.3, 2.8, 1.0, 0.6, size=10)

# TEMPO Downlooking Dual Vector Indicators
arrow_sat1 = slide.shapes.add_shape(MSO_SHAPE.DOWN_ARROW, Inches(5.4), Inches(1.7), Inches(0.2), Inches(0.6))
arrow_sat2 = slide.shapes.add_shape(MSO_SHAPE.DOWN_ARROW, Inches(5.7), Inches(1.7), Inches(0.2), Inches(0.6))
for arr in [arrow_sat1, arrow_sat2]:
    arr.fill.solid(); arr.fill.fore_color.rgb = ARROW_BLUE; arr.line.fill.background()

# Callout Bubble Indicators
bubble_sun = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.6), Inches(1.5), Inches(1.6), Inches(1.0))
bubble_sun.fill.solid(); bubble_sun.fill.fore_color.rgb = WHITE
add_clean_label("P-3\nSounding\nRange", 0.7, 1.6, 1.4, 0.8, size=10, center=True)

panel_meas = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(5.2), Inches(2.4), Inches(2.0), Inches(0.9))
panel_meas.fill.solid(); panel_meas.fill.fore_color.rgb = WHITE
add_clean_label("total column\nsunphotometer", 5.25, 2.45, 1.9, 0.8, size=10, center=True)

# Data Assimilation Cyclic Framework Block
loop_arrow = slide.shapes.add_shape(MSO_SHAPE.BLOCK_ARC, Inches(9.2), Inches(1.5), Inches(3.5), Inches(3.5))
loop_arrow.fill.solid(); loop_arrow.fill.fore_color.rgb = ARROW_BLUE; loop_arrow.line.fill.background()

add_clean_label("WRF-Chem\nCTM", 9.9, 2.8, 2.2, 0.8, size=16, bold=True, center=True)
add_clean_label("Model-Measurement\nSynthesis", 9.7, 3.8, 2.6, 0.5, size=10, bold=False, center=True)

# --- 6. Cartesian Profiling Graphs Axis Dashboards ---
def draw_profile_axis(x_start, name, axis_max_lbl):
    frame = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x_start), Inches(1.2), Inches(1.2), Inches(5.3))
    frame.fill.solid(); frame.fill.fore_color.rgb = WHITE; frame.line.color.rgb = TEXT_DARK
    
    add_clean_label(name, x_start, 0.6, 1.2, 0.4, size=11, bold=True, center=True)
    add_clean_label("Profile", x_start, 0.85, 1.2, 0.3, size=10, bold=True, center=True)
    
    # Altitude markings along left canvas boundary of individual chart
    add_clean_label("10", x_start - 0.25, 2.2, 0.3, 0.3, size=9, bold=False)
    add_clean_label("5", x_start - 0.25, 4.2, 0.3, 0.3, size=9, bold=False)
    add_clean_label("2", x_start - 0.25, 5.4, 0.3, 0.3, size=9, bold=False)
    add_clean_label("0", x_start - 0.25, 6.3, 0.3, 0.3, size=9, bold=False)
    add_clean_label(axis_max_lbl, x_start + 0.8, 6.5, 0.4, 0.3, size=8, bold=False)

# Render Profiles Side-by-Side on Right Dashboard Area
draw_profile_axis(8.8, "NO2", "4")
draw_profile_axis(10.3, "Aerosol", "1")
draw_profile_axis(11.8, "O3", "80")

# Save complete editable graphic layout presentation
prs.save("Atmospheric_Assimilations_Complete_Vector.pptx")
print("Full interactive diagram vectors generated cleanly as 'Atmospheric_Assimilations_Complete_Vector.pptx'.")


# In[ ]:




