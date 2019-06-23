# Analyzing the Effects of Space Flight on Telomere Length Dynamics in NASA Astronauts
  
In 2024 NASA will execute Project Artemis, sending humans back to the moon to establish a permanent lunar base. However, space is inherently dangerous to human health. As part of my Ph.D, I'm researching how spaceflight impacts human health, and whether these impacts could potentially comprise current (Project Artemis) or future missions, ala Mars and beyond. Specifically, I'm examining how time aboard the International Space Station affects telomeres, the ends of human DNA, and the stability of DNA, for NASA's astronauts. My research takes the *first look at the changes to telomeres in unrelated astronauts as a result of spaceflight*, informing NASA policy and approach for current and future missions.

This repository details the context of my research in the wider pursuit of space exploration; a brief overview of my laboratory methods; and an extensive walkthrough of my analysis using Python on telomere length data in astronauts. Why Python? The first pass of my analysis on the astronaut dataset was done in Excel and took ~1 month; my Python script takes about 30 seconds. This repository may serve as an accessory to our upcoming publications. The walkthrough below is intended for those interested in the science and in the code. A Jupyter notebook will be provided as an accessory. 

Please feel free to contact me.

**Contact:**  
Jared Luxton  
jLuxton@colostate.edu

<p align="center">
<a href="url">
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Python-logo-notext.svg/200px-Python-logo-notext.svg.png" height="150"> 
<img src="https://cdn1.medicalnewstoday.com/content/images/articles/319/319971/space-explorer.jpg" height="150">
<img src="https://abm-website-assets.s3.amazonaws.com/rdmag.com/s3fs-public/embedded_image/2017/04/telomere-chromosome-stock.jpg" height="150">
</a>
</p>
&nbsp;
&nbsp;   

## Table of Contents:
* [Background: Health Risks and Obstacles to Space Flight](#background) 
* [Approach: Identifying Risks by Measuring Telomere Lengths](#approach)
* [Methods: Blood Collections, Cell Culture, Telomere Measurements](#methods)
* [Data Cleaning: Handling Telomere Length Data](#data-cleaning)
* [Data Analysis: Visualization and Statistics](#data-analysis)
* [Conclusions](#conclusions)

&nbsp;    

## Background 
**Health Risks and Obstacles to Space Flight**\
Did you know that NASA is sending humans to the Moon in 2024? Yes! And not only that: this mission is the first of many that will develop permanent lunar colonies and provide a bridge to exploring Mars and beyond. Dubbed Artemis, this NASA project entails sending the *first* woman (and another man) to the lunar surface and the development of a *permanent lunar outpost called the Gateway* orbiting the Moon. The objectives undertaken by Artemis are part of NASA's overarching goal (and humanity's common dream) for humans to explore our solar system; Mars and beyond. 

The immediate challenges facing Artemis are substantial in terms of technology and health considerations for the astronauts. Even as we approach Project Artemis in 2024, the short- and long-term health effects of spaceflight, especially those from chronic exposure to galactic cosmic rays, a type of radiation unique to space and not found on Earth, remain relatively unknown. 

Galactic cosmic rays (GCRs) are highly energetic particles hurtling through space at nearly the speed of light. Though a rare event, when GCRs strike human cells they shred all cellular contents in their path, including DNA. This damage accumulates over time, and could lead to degeneration of tissues and cancer. Currently, we simply don't understand how much cellular damage humans accumulate in space, and how much it increases cancer risk. Not understanding these issues makes addressing them impossible. My research directly addresses these issues by examing how spaceflight effects telomeres (the ends of DNA) and DNA stability for NASA astronauts aboard the International Space Station. 

&nbsp; 

## Approach 
**Identifying Risks by Measuring Telomere Lengths**\
Telomeres are repetitive sequences of DNA covered by protein found at the very ends of DNA. Telomeres shorten with each cell division and thus shorten as we age. When the telomeres in a cell reach a critically short length, the cell will die or persist in a state which damages neighboring cells (termed senescence). Cell death resulting from telomeres shortening too quickly will lead to age-related diseases, i.e cancer.  Environmental exposures, including space radiation, air pollution, stress, inflammation, and others can all contribute to telomere shortening and thus age-related diseases - cancer. Telomeres therefore link environmental exposures with age-related diseases. By measuring telomere length over a period of time which involves environmental stressors and exposures, the telomere length changes can be used to quantify the short- and long-term effects of that experience in terms of cancer risk and disease. This is what we've done with the astronauts.

&nbsp; 

## Methods
**Blood Collections, Cell Culture, Telomere Measurements**\
We monitored the telomere lengths in 11 unrelated astronauts at pre-, mid-, and post-spaceflight timepoints aboard the International Space Station (where available). In all, we have telomere data for about seven timepoint samples for each astronaut. For our analyses, we directly monitored the lengths of *all individual telomeres* in each cell for each timepoint for each patient; we also have the telomere length means for those timepoints.

To measure telomere length in astronauts, we used a noninvasive approach for sample collection and analysis. Blood was taken from astronauts at pre-, mid- (yes, blood was drawn aboard the International Space Station, sent down on the Soyuz capsule to Texas, and mailed to us), and post-spaceflight timepoints. From this blood we specifically cultured white blood cells (ala 'T-cells'). By culturing and using only white blood cells for our telomere length measurements, we ensured the homogenity of our measurements. 

&nbsp; 

## Data Cleaning 
**Handling Telomere Length Data**\
Here's what you've been reading for. *The Python good stuff.* 

```python
s = "Python syntax highlighting"
print('HELLO WORLD!')
'''

&nbsp; 

## Data Analysis
**Visualization and Statistics**\
...

&nbsp; 

## Conclusions
**Highlights and Final Thoughts**\
...

&nbsp; 
&nbsp; 
