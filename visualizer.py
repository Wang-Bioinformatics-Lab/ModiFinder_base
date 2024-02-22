import xlsxwriter
import os
from . import utils_n as utils
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
import cairosvg
import pandas as pd
import urllib.parse
import matplotlib.pyplot as plt
import uuid
import math
import numpy as np
import matplotlib.image as mpimg
import io

def make_url(base_url = "http://localhost:8050/", USI1=None, USI2=None, SMILES1=None, SMILES2=None, args=None):
    query_params = {k: v for k, v in {
        "USI1": USI1,
        "USI2": USI2,
        "SMILES1": SMILES1,
        "SMILES2": SMILES2
    }.items() if v is not None}
    
    if args is not None:
        query_params.update(args)

    url = base_url + "?" + urllib.parse.urlencode(query_params)
    return url

def table_to_xlsx(data, path):
    """
    Converts a dataframe to an Excel file.
    """
    wb = xlsxwriter.Workbook(path)
    my_format = wb.add_format()
    my_format.set_align('vcenter')
    ws = wb.add_worksheet('Sheet 1')
    col = 0
    temp_file_names = []
    for key in data:
        if 'image' not in key and key != 'structure':
            # print(key, data[key])
            ws.write(0, col, key)
            ws.write(1, col, data[key])
            ws.set_column(col, col, None, my_format)
            col += 1
        else:
            ws.write(0, col, key)
            # write svg to file
            filename = str(uuid.uuid4())
            with open("{}.svg".format(filename), "w") as f:
                f.write(data[key])
                temp_file_names.append("{}.svg".format(filename))
            # svg image to bitmap
            cairosvg.svg2png(url="{}.svg".format(filename), write_to="{}.png".format(filename), parent_width=250, parent_height=200)
            temp_file_names.append("{}.png".format(filename))

            # make image fit in cell
            ws.set_column_pixels(col, col, 150)
            ws.set_row_pixels(1, 150)
            ws.insert_image(1, col, "{}.png".format(filename), {'x_scale': 0.12, 'y_scale': 0.12, 'x_offset': 1, 'y_offset': 1})
            col += 1
    
    wb.close()
    for filename in temp_file_names:
        try:
            os.remove(filename)
        except:
            pass

def table_to_tsv(data, path_to_tsv, path_to_img_folder):
    import uuid

    unique_key = str(uuid.uuid4())

    path_to_svg = os.path.join(path_to_img_folder, unique_key + ".svg")

    data["img_path"] = path_to_svg

    df = pd.DataFrame([data])
    df = df[['# matched peaks', '# unshifted peaks', '# shifted peaks', 'usi1', 'usi2', 'structure1', 'img_path']]
    
    # TODO: Adding structures and USIs

    df.to_csv(path_to_tsv, sep='\t', index=False)

    # Writing out the images to file
    with open(path_to_svg, "w") as f:
        f.write(data["image"])


def molToSVG(mol, substructure=None, highlightModificationSites=False):
    """
    Converts a molecule to an SVG image.
    """
    if substructure is not None:
        if (type(substructure) == str):
            substructure = Chem.MolFromSmarts(substructure)
        if not mol.HasSubstructMatch(substructure):
            if substructure.HasSubstructMatch(mol):
                return molToSVG(substructure, mol, highlightModificationSites)
            else:
                print("there")
                d2d = Chem.Draw.MolDraw2DSVG(1250,1200)
                d2d.DrawMolecule(mol)
                d2d.FinishDrawing()
                return d2d.GetDrawingText()
        
        hitAtoms, hitBonds = utils.getHitAtomsAndBonds(mol, substructure)
        highlightAtoms = []
        for atom in mol.GetAtoms():
            if atom.GetIdx() not in hitAtoms[0]:
                highlightAtoms.append(atom.GetIdx())
        
        highlightbonds = []
        for bond in mol.GetBonds():
            # if both side are in highlightAtoms then add to highlightbonds
            if bond.GetBeginAtomIdx() in highlightAtoms or bond.GetEndAtomIdx() in highlightAtoms:
                highlightbonds.append(bond.GetIdx())
        

        if highlightModificationSites:
            # print(mol.GetNumAtoms(), substructure.GetNumAtoms(), "To check the sizes!")
            modifications = utils.calculateModificationSites(mol, substructure)
            # print(modifications)
            d2d = Chem.Draw.MolDraw2DSVG(1250,1200)
            colors = dict()
            for hitAtom in highlightAtoms:
                colors[hitAtom] = (1,0.2,0.2)
            for atom in mol.GetAtoms():
                if atom.GetIdx() in modifications:
                    colors[atom.GetIdx()] = (0.2,0.2,1)
            d2d.DrawMolecule(mol,highlightAtoms=highlightAtoms, highlightBonds=highlightbonds, highlightAtomColors=colors)
        else:
            d2d = Chem.Draw.MolDraw2DSVG(1250,1200)
            d2d.DrawMolecule(mol,highlightAtoms=highlightAtoms, highlightBonds=highlightbonds)
    else:
        d2d = Chem.Draw.MolDraw2DSVG(1250,1200)
        d2d.DrawMolecule(mol)
    d2d.FinishDrawing()



    return d2d.GetDrawingText()

def draw_png(mol):
    d2d = Draw.MolDraw2DCairo(1200, 1200)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    png = d2d.GetDrawingText()
    img = mpimg.imread(io.BytesIO(png))
    return img

def highlightScores(mol, scores, add_labels = False, shrink_labels = False, label_size = 0.5):
    if min(scores) < 0:
        scores = [x + abs(min(scores)) for x in scores]
    if max(scores) > 0:
        vals = [x/max(scores) for x in scores]
    else:
        vals = [0 for x in scores]
    d2d = Draw.MolDraw2DSVG(1250,1200)
    colors = dict()
    for i in range(0, mol.GetNumAtoms()):
        # heat map coloring
        if vals[i] == 0:
            colors[i] = (vals[i], 0.2*(1-vals[i]), 1-vals[i], 0.4)
        else:
            colors[i] = (vals[i], 0.2*(1-vals[i]), 1-vals[i], vals[i]*0.3 + 0.65)
    if add_labels:
        d2d.drawOptions().annotationFontScale = label_size
        for atom in mol.GetAtoms():
            lbl = str(round(scores[atom.GetIdx()], 2))
            if shrink_labels:
                if scores[atom.GetIdx()] == 0:
                    lbl = ""
                else:
                    lbl = str(int(round(scores[atom.GetIdx()], 2)*100))
            atom.SetProp('atomNote',lbl)
            # set font size for labels
    d2d.DrawMolecule(mol, highlightAtoms=list(range(mol.GetNumAtoms())), highlightAtomColors=colors, highlightBonds=[])
    d2d.FinishDrawing()
    def draw_gradient_svg(length, diam, ax = 0, steps = 100, fontSize = 10):
        '''
        ax = 0: x axis
        ax = 1: y axis
        '''

        width = length
        height = diam
        if ax == 1:
            width = diam
            height = length
        
        svg = """<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg">""".format(width, 1200)
        rectWidth = width/steps
        rectHeight = height
        if ax == 1:
            rectWidth = width
            rectHeight = height/steps


        # else:
        #     # add text to svg rotated
        #     svg += """<text x="{}" y="{}" fill="black" transform="rotate(90, {}, {})">{}</text>""".format(0, height/2, 0, height/2, "low likelihood")
        #     svg += """<text x="{}" y="{}" fill="black" transform="rotate(90, {}, {})">{}</text>""".format(0, width-100, 0, width-100, "high likelihood")

        for i in range(steps, -1, -1):
            if ax == 0:
                svg += """<rect x="{}" y="{}" width="{}" height="{}" fill="rgba({}, {}, {}, {})"/>""".format(i*rectWidth, 1200-diam, rectWidth*1.01, rectHeight, (i/steps)*255, 0.2*(1-i/steps)*255,(1-i/steps)*255,  ((i/steps)*0.3) + 0.65)
            else:
                svg += """<rect x="{}" y="{}" width="{}" height="{}" fill="rgba({}, {}, {}, {})"/>""".format(0, i*rectHeight, rectWidth(1.01), rectHeight, (i/steps)*255, 0.2*(1-i/steps)*255, (1-i/steps)*255,  ((i/steps)*0.3) + 0.65)
        
        if ax == 0:
            svg += """<text x="{}" y="{}" font-size="{}px" fill="white">{}</text>""".format(5, 1200-diam + (height+fontSize/2)/2, fontSize, "low likelihood")
            # get width of text
            text = "high likelihood"
            svg += """<text x="{}" y="{}" font-size="{}px" fill="white">{}</text>""".format(width-len(text)*fontSize/2.4, 1200-diam + (height+fontSize/2)/2, fontSize, text)
        else:
            text = "high likelihood"
            # change font size
            
            svg += """<text x="{}" y="{}" font-size="{}px" fill="white" transform="rotate(90, 0, 0)">{}</text>""".format(5, -(width-fontSize)/2, fontSize, "low likelihood")
            svg += """<text x="{}" y="{}" font-size="{}px" fill="white" transform="rotate(90, 0, 0)">{}</text>""".format(height-len(text)*fontSize/1.5, -(width-fontSize)/2, fontSize, text)

        svg += "</svg>"
        return svg
    
    svgText = d2d.GetDrawingText()
    svgGradient = draw_gradient_svg(1250, 50, ax = 0, steps = 30, fontSize = 40)

    # concatenate two svg files
    svgText = svgText.replace("</svg>", svgGradient + "</svg>")

    if add_labels:
        # remove labels
        for atom in mol.GetAtoms():
            # lbl = '%.2f'%(scores[atom.GetIdx()])
            atom.SetProp('atomNote', '')
    
    return svgText

def highlightMolIndices(mol, hitAtoms, hitBonds = None, resolution = None):
    if resolution is not None:
        d2d = Draw.MolDraw2DSVG(resolution[0], resolution[1])
    else:
        d2d = Draw.MolDraw2DSVG(1250,1200)
    if hitBonds is None:
        hitBonds = []
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in hitAtoms and bond.GetEndAtomIdx() in hitAtoms:
                hitBonds.append(bond.GetIdx())
    d2d.DrawMolecule(mol,highlightAtoms=hitAtoms, highlightBonds=hitBonds)
    d2d.FinishDrawing()
    svgText =  d2d.GetDrawingText()
    return svgText

def draw_alignment(peaks1, peaks2, matched_peaks, shift = 0.1, show_text = False, show_lines = True, scale = 1, ax = None, flip = True, x_lim = None, legend_loc = 'upper right', text_size = 1):

    if ax is None:
        fig, ax = plt.subplots(figsize=(20*scale, 5*scale))

    shifted, unshifted = utils.separateShifted(matched_peaks, peaks1, peaks2)
    max_y = max([_[1] for _ in peaks1])
    peaks1 = [(peaks1[i][0], peaks1[i][1]/max_y) for i in range(len(peaks1))]

    max_y = max([_[1] for _ in peaks2])
    peaks2 = [(peaks2[i][0], peaks2[i][1]/max_y) for i in range(len(peaks2))]

    max_x = max([_[0] for _ in peaks1])
    max_x = max([max_x, max([_[0] for _ in peaks2])])

    ax.set_xlim(0, max_x+5)

    print('debugging', peaks2)
    
    molShifted = [i[0] for i in shifted]
    molUnshifted = [i[0] for i in unshifted]
    modifiedShifted = [i[1] for i in shifted]
    modifiedUnshifted = [i[1] for i in unshifted]

    # draw peaks in site locator main compound
    # make shifted peaks red
    # make unshifted peaks blue

    x = [_[0] for _ in peaks1 ]
    y = [ _[1] for _ in peaks1 ]
    flip_shift = 0
    if not flip:
        flip_shift = 1
    shift_array = [shift + flip_shift for _ in peaks1 ]
    ax.bar(x, shift_array, width=0.3*scale, color="white", alpha = 0)
    ax.bar(x, y, width=0.3*scale, color="gray", bottom=shift_array)

    x_shifted = [peaks1[_][0] for _ in molShifted]
    y_shifted = [ peaks1[_][1] for _ in molShifted]
    shift_array = [shift + flip_shift for _ in molShifted]
    ax.bar(x_shifted, shift_array, width=0.3*scale, color="white", alpha = 0)
    ax.bar(x_shifted, y_shifted, width=0.3*scale, color="red", bottom=shift_array)

    x_unshifted = [peaks1[_][0] for _ in molUnshifted]
    y_unshifted = [ peaks1[_][1] for _ in molUnshifted]
    shift_array = [shift + flip_shift for _ in molUnshifted]
    ax.bar(x_unshifted, shift_array, width=0.3*scale, color="white", alpha = 0)
    ax.bar(x_unshifted, y_unshifted, width=0.3*scale, color="blue", bottom=shift_array)

    #plot modified peaks as reversed
    x = [_[0] for _ in peaks2]
    if not flip:
        y = [_[1] for _ in peaks2]
        y_shifted = [peaks2[_][1] for _ in modifiedShifted]
        y_unshifted = [peaks2[_][1] for _ in modifiedUnshifted]
    else:
        y = [-_[1] for _ in peaks2]
        y_shifted = [-peaks2[_][1] for _ in modifiedShifted]
        y_unshifted = [-peaks2[_][1] for _ in modifiedUnshifted]

    ax.bar(x, y, width=0.3*scale, color="gray")
    x_shifted = [peaks2[_][0] for _ in modifiedShifted]
    ax.bar(x_shifted, y_shifted, width=0.3*scale, color="red")
    x_unshifted = [peaks2[_][0] for _ in modifiedUnshifted]
    ax.bar(x_unshifted, y_unshifted, width=0.3*scale, color="blue")

    if show_lines:
        for peak in shifted:
            x1 = peaks1[peak[0]][0]
            x2 = peaks2[peak[1]][0]

            if x_lim is not None:
                if x1 < x_lim[0] or x1 > x_lim[1] or x2 < x_lim[0] or x2 > x_lim[1]:
                    continue

            y1 = shift + flip_shift
            y2 = 0
            # draw line with text in the middle
            ax.plot([x1, x2], [y1, y2], color="red", linewidth=0.5*scale, linestyle="--")
            val = abs(peaks1[peak[0]][0] - peaks2[peak[1]][0])
            val = round(val, 2)
            if show_text:
                ax.text((x1 + x2) / 2, shift/2, str(val), fontsize=10, horizontalalignment='center')

        for peak in unshifted:
            x1 = peaks1[peak[0]][0]
            x2 = peaks2[peak[1]][0]
            y1 = shift + flip_shift
            y2 = 0
            ax.plot([x1, x2], [y1, y2], color="blue", linewidth=0.5*scale, linestyle="--")

    #draw horizontal line
    ax.plot([2, max_x], [shift + flip_shift, shift + flip_shift], color="gray", linewidth=0.3*scale, linestyle="-")
    ax.plot([2, max_x], [0, 0], color="gray", linewidth=0.2*scale, linestyle="-")

    # custom y axis ticks
    y_ticks1 = [i/10 + shift + flip_shift for i in range(0, 11, 2)]
    if flip:
        y_ticks2 = [-i/10 for i in range(0, 11, 2)]
    else:
        y_ticks2 = [i/10 for i in range(0, 11, 2)]
    # reverse y ticks2
    y_ticks2 = y_ticks2[::-1]
    y_ticks =  y_ticks2 + y_ticks1 
    y_ticks_labels1 = [i for i in range(0, 110, 20)]
    y_ticks_labels2 = [i for i in range(0, 110, 20)]
    # reverse y ticks2
    y_ticks_labels2 = y_ticks_labels2[::-1]
    y_ticks_labels = y_ticks_labels2 + y_ticks_labels1
    ax.set_yticks(y_ticks, y_ticks_labels)

    # set font size
    ax.tick_params(axis='both', which='major', labelsize=20*text_size)

    # add legend
    ax.plot([], [], color="red", linewidth=2*scale, linestyle="-", label="Shifted Matched Peaks")
    ax.plot([], [], color="blue", linewidth=2*scale, linestyle="-", label="Unshifted Matched Peaks")
    ax.plot([], [], color="gray", linewidth=2*scale, linestyle="-", label="Unmatched Peaks")
    
    # legend font size and position
    ax.legend(prop={'size': 20*text_size}, loc=legend_loc)

    # set x range
    if x_lim is not None:
        ax.set_xlim(x_lim[0], x_lim[1])

    # remove top and right borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return ax


def render_svg_to_png(svg_content):
    """Render SVG content to PNG using cairosvg."""
    s_png = cairosvg.svg2png(bytestring=svg_content)
    s_img = mpimg.imread(io.BytesIO(s_png))
    return s_img

def draw_frags_of_peak(modSiteLocator, peak_index, fig = None):
    structures, locations, frags = modSiteLocator.get_structures_by_peak_index(peak_index)
    subplotCols = 2
    subplotsRows = math.ceil(len(locations)/subplotCols)
    if fig is None:
        fig, ax = plt.subplots(subplotsRows, subplotCols, figsize=(100, 100))
    else:
        ax = []
        for i in range(subplotsRows):
            row = []
            for j in range(3):
                row.append(fig.add_subplot(subplotsRows, subplotCols, i*subplotCols+j+1))
            ax.append(row)
        ax = np.array(ax)
        if subplotsRows == 1:
            ax = ax[0]
    for i, location in enumerate(locations):
        svgText = highlightMolIndices(modSiteLocator.main_compound.structure, location)
        img = render_svg_to_png(svgText)
        if subplotsRows == 1:
            ax[i].imshow(img)
            # set title to fragment
            ax[i].set_title("fragment " + str(frags[i]))
            ax[i].axis('off')
        else:
            ax[i//subplotCols, i%subplotCols].imshow(img)
            # set title to fragment
            ax[i//subplotCols, i%subplotCols].set_title("fragment " + str(frags[i]))
        
            # remove axis
            ax[i//subplotCols, i%subplotCols].axis('off')
    
    # make the ax tight
    fig.tight_layout()
