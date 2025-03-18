# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:20:18 2025

@author: jhodges
"""

import glob, os

def checkCaption(caption):
    
    # Punctuation not allowed at the end of the short name
    captionTocDisallowedPunctuation = '.!?'
    
    # Missing short name handling
    missing_toc_name_style = 1 # 0 - do not flag, 1 - flag as warning, 2 - flag as error
    
    # Initialize outtxt
    outtxt = ''
    
    # Check for presence of short name, needed for check for citation in TOC regardless
    if '[' not in caption.split('{')[0]: 
        notoc_name = True
        short_name = ' '
    else:
        notoc_name = False
        short_name = caption.split('[')[1].split(']')[0]
    
    if notoc_name and missing_toc_name_style == 1:
        outtxt = outtxt + "WARNING, %s caption, %s does not have a TOC name"%(file,caption) + "\n"
    elif notoc_name and missing_toc_name_style == 2:
        outtxt = outtxt + "ERROR, %s caption, %s does not have a TOC name"%(file,caption) + "\n"
    
    # Check if short name ends in disallowed puncuation
    if short_name[-1] in captionTocDisallowedPunctuation:
        outtxt = outtxt + "WARNING, %s caption, %s TOC name ends in '%s'"%(file,caption, short_name[-1]) + "\n"
    
    # Check if citation is included in the name used in TOC
    if notoc_name and '\\cite' in caption:
        outtxt = outtxt + "ERROR, %s citation in caption, %s"%(file,caption) + "\n"
    elif '\\cite' in short_name:
        outtxt = outtxt + "ERROR, %s citation in caption, %s"%(file,caption) + "\n"
        
    return outtxt

def check_disallowed_commands(txt, file):
    disallowed_commands = ['\\bf{','\\tt{']
    outtxt = ''
    for cmd in disallowed_commands:
        split = txt.split(cmd)
        if len(split) > 1:
            for j in range(1, len(split)):
                line_count = len(split[j-1].split('\n'))+1
                outtxt = outtxt + "ERROR, %s %s located at line %d\n"%(file, cmd, line_count)
    return outtxt

texfiles = glob.glob('..'+os.sep+'*'+os.sep+'*.tex')

outtxt = '\n'
for i in range(0, len(texfiles)):
    file = texfiles[i]
    
    with open(file, 'r') as f:
        txt = f.read()
    
    # Check figures
    figs = txt.split('begin{figure}')
    for j in range(1, len(figs)):
        fig = figs[j].split('end{figure}')[0]
        captions = fig.split('\\caption')
        for k in range(1, len(captions)):
            caption = captions[k].split('}')[0] + '}'
            outtxt = outtxt + checkCaption(caption)
            
    # Check tables
    tabs = txt.split('begin{table}')
    for j in range(1, len(tabs)):
        tab = tabs[j].split('end{table}')[0]
        captions = tab.split('\\caption')
        for k in range(1, len(captions)):
            caption = captions[k].split('}')[0]
            outtxt = outtxt + checkCaption(caption)
    
    # Check disallowed commands
    outtxt = outtxt + check_disallowed_commands(txt, file)

if len(outtxt) > 1:
    print("Warnings identified in the manual check:")
    print(outtxt)

with open('check_output.txt', 'w') as f:
    f.write(outtxt)