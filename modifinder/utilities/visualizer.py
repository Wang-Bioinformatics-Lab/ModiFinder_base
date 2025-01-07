"""
GNPS Utils - Visualizer Module

This module provides functionality to visualize data from GNPS.

Author: Shahneh
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdFMCS
from rdkit.Chem import Draw
from modifinder.convert import to_spectrum
from modifinder.utilities.network import *
from modifinder.utilities.gnps_types import *
from modifinder.utilities.general_utils import *
import modifinder.utilities.mol_utils as mu
# import modifinder.utilities.spectra_utils as su
# from modifinder.alignment import _cosine_fast
import io
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from PIL import ImageDraw, ImageFont, Image
from matplotlib.patches import ConnectionPatch
import os
from io import BytesIO

def draw_molecule(mol, output_type='png', font_size = None, label=None, label_font_size = 20, label_color = (0,0,0), label_position = 'top', **kwargs):
    """
    Draw a molecule using RDKit
    
    Parameters
    ----------
    mol : rdkit molecule or str
        rdkit molecule or str for SMILES or InChI or GNPS identifier (USI or Accession)
    output_type : str
        type of output (png or svg)
    font_size : int, optional (default=None)
        font size for the labels
    label : str, optional (default=None)
        label for the molecule
    label_font_size : int, optional (default=20)
        font size for the label
    label_color : tuple, optional (default=(0,0,0))
        color of the label
    label_position : str, optional (default='top')
        position of the label (top or bottom)
    kwargs : dict
        additional arguments for drawing the molecule in rdkit
        like highlightAtoms, highlightAtomColors, highlightBonds, highlightBondColors, highlightAtomRadii, etc.
    
    Returns
    -------
    img : numpy array or str
        image of the molecule
    
    Examples
    -------
    .. code-block:: python
    
        import modifinder.utilities.visualizer as mf_vis
        from matplotlib import pyplot as plt

        img = mf_vis.draw_molecule('CN1C=NC2=C1C(=O)N(C(=O)N2C)C', output_type='png', label="Caffeine")
        plt.imshow(img)
        plt.axis('off')
        plt.show()
    
    .. image:: ../_static/draw_molecule1.png
        :width: 300px
    
    .. code-block:: python
    
        import modifinder.utilities.visualizer as mf_vis
        from matplotlib import pyplot as plt
        from rdkit import Chem

        def mol_with_atom_index(mol):
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx())
            return mol
        mol = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
        mol = mol_with_atom_index(mol)
        highlightAtoms = {1, 3, 11, 8}
        img = mf_vis.draw_molecule(mol, output_type='png', label="Caffeine", highlightAtoms=highlightAtoms)
        plt.imshow(img)
        plt.axis('off')
        plt.show()
    
    .. image:: ../_static/draw_molecule2.png
        :width: 300px
    """
    # TODO: test svg

    molecule = mu._get_molecule(mol)
    if molecule is None:
        raise ValueError("Molecule not found")
    
    if type(output_type) != str:
        raise ValueError("Output type should be a string")
    output_type = output_type.lower()
    if output_type not in ["png", "svg"]:
        raise ValueError("Output type should be either png or svg")
    
    # if size is not provided, use default size
    x_dim, y_dim = kwargs.get("size", (400, 400))
    x_dim = kwargs.get("x_dim", x_dim)
    y_dim = kwargs.get("y_dim", y_dim)

    if font_size is None:
        font_size = x_dim // 20

    extra_info = kwargs.get("extra_info", {})

    draw_kwargs = {}
    for key in ['highlightAtoms', 'highlightAtomColors', 'highlightBonds', 'highlightBondColors', 'highlightAtomRadii']:
        if key in kwargs:
            draw_kwargs[key] = kwargs[key]

    if output_type == "png":
        d2d = Draw.MolDraw2DCairo(x_dim, y_dim)
        d2d.drawOptions().minFontSize = font_size
        d2d.drawOptions().maxFontSize = font_size
        if "annotation_scale" in extra_info:
            d2d.drawOptions().annotationFontScale = extra_info["annotation_scale"]        

        d2d.DrawMolecule(molecule, **draw_kwargs)
        d2d.FinishDrawing()
        png = d2d.GetDrawingText()
        img = mpimg.imread(io.BytesIO(png))

        if label:
            img = Image.fromarray((img*255).astype(np.uint8))
            draw = ImageDraw.Draw(img)
            path_of_this_file = os.path.abspath(os.path.dirname(__file__))
            font = ImageFont.truetype(os.path.join(path_of_this_file, "fonts/NotoSans-Regular.ttf"), label_font_size)
            if label_position == 'top':
                draw.text((x_dim//2, 0), label, label_color, font=font, anchor="mt")
            else:
                draw.text((x_dim//2, y_dim), label, label_color, font=font, anchor="mb")
            img = np.array(img)

        if "show_legend" in extra_info and extra_info["show_legend"]:
            img = _add_heatmap_legend(img, extra_info["scores"], kwargs["highlightAtomColors"],
                               x_dim, y_dim, extra_info["legend_width"], extra_info["legend_font"], output_type)
        return img
    else:
        d2d = Chem.Draw.MolDraw2DSVG(x_dim, y_dim)
        d2d.SetFontSize(font_size)
        if "annotation_scale" in extra_info:
            d2d.drawOptions().annotationFontScale = extra_info["annotation_scale"]
        d2d.DrawMolecule(molecule, **draw_kwargs)
        d2d.FinishDrawing()
        svg = d2d.GetDrawingText()
        if "show_legend" in extra_info and extra_info["show_legend"]:
            svg = _add_heatmap_legend(svg, extra_info["scores"], kwargs["highlightAtomColors"],
                              x_dim, y_dim, extra_info["legend_width"], extra_info["legend_font"], output_type)
        return svg


def draw_modifications(mol1, mol2, output_type='png', show_legend = True, legend_font = 15, legend_position = None, highlight_common = True, highlight_added = True, highlight_removed=True, modification_only = False, **kwargs):
    """
    Draw the modifications from molecule 1 to molecule 2
    
    Parameters
    ----------
    mol1 : rdkit molecule, str
        rdkit molecule or str for SMILES or InChI or GNPS identifier (USI or Accession)
    mol2 : rdkit molecule
        rdkit molecule or str for SMILES or InChI or GNPS identifier (USI or Accession)
    output_type : str, optional (default='png')
        type of output (png or svg)
    show_legend : bool, optional (default=True)
        show the legend or not
    legend_font : int, optional (default=15)
        font size of the legend
    legend_position : tuple, optional (default=None)
        position of the legend
    highlight_common : bool, optional (default=True)
        highlight the common atoms
    highlight_added : bool, optional (default=True)
        highlight the added atoms
    highlight_removed : bool, optional (default=True)
        highlight the removed atoms
    modification_only : bool, optional (default=False)
        only highlight the modification edges, if highlight_removed or highlight_added is False, this will only show the enabled ones
    
    Returns
    -------
    img : numpy array or str
        image of the modification
    
    Examples
    -------
    .. code-block:: python
    
        import modifinder.utilities.visualizer as mf_vis
        from matplotlib import pyplot as plt
        from rdkit import Chem

        smiles1, smiles2 = 'N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(O)=O)C(O)=O', 'CCCCCCSCC(CNCC(=O)O)NC(=O)CCC(C(=O)O)N'
        mol1 = mf_vis.draw_molecule(smiles1, label="mol1")
        mol2 = mf_vis.draw_molecule(smiles2, label="mol2")
        modification = mf_vis.draw_modifications(smiles1, smiles2, label="modifications")
        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        ax[0].imshow(mol1)
        ax[1].imshow(mol2)
        ax[2].imshow(modification)
        for a in ax:
            a.axis('off')
        plt.show()
    
    .. image:: ../_static/draw_modifications.png
    """
    mol1, mol2 = mu._get_molecules(mol1, mol2)
    result = mu.get_transition(mol1, mol2)

    highlight_color_removed = (0.8, 0.1, 0.1, 0.2)
    highlight_color_removed_dark = (1, 0.1, 0, 0.8)
    highlight_color_added = (0.1, 0.8, 0.1, 0.2)
    highlight_color_added_dark = (0, 0.7, 0.2, 0.7)
    highlight_color_common = (0.1, 0.1, 0.8, 0.2)

    highlight_atoms = []
    highlight_bonds = []
    highlight_atoms_colors = dict()
    highlight_bonds_colors = dict()

    if modification_only:
        highlight_common = False

    if highlight_common:
        for atom in result['common_atoms']:
            highlight_atoms.append(atom)
            highlight_atoms_colors[atom] = highlight_color_common
    
    if highlight_added:
        for atom in result['added_atoms']:
            highlight_atoms.append(atom)
            highlight_atoms_colors[atom] = highlight_color_added

    if highlight_removed:
        for atom in result['removed_atoms']:
            highlight_atoms.append(atom)
            highlight_atoms_colors[atom] = highlight_color_removed

    modification_added_edges = result['modified_added_edges_bridge'] + result['modified_added_edges_inside']
    modification_removed_edges = result['modified_removed_edges_bridge'] + result['modified_removed_edges_inside']

    for bondIdx in range(result['merged_mol'].GetNumBonds()):
        bond = result['merged_mol'].GetBondWithIdx(bondIdx)
        pair1 = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        pair2 = (bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())
        if highlight_common:
            if pair1 in result['common_bonds'] or pair2 in result['common_bonds']:
                highlight_bonds_colors[bondIdx] = highlight_color_common
                highlight_bonds.append(bondIdx)
        if highlight_added:
            if pair1 in modification_added_edges or pair2 in modification_added_edges:
                highlight_bonds_colors[bondIdx] = highlight_color_added_dark
                highlight_bonds.append(bondIdx)
            elif (not modification_only) and (pair1 in result['added_edges'] or pair2 in result['added_edges']):
                highlight_bonds_colors[bondIdx] = highlight_color_added
                highlight_bonds.append(bondIdx)
        if highlight_removed:
            if pair1 in modification_removed_edges or pair2 in modification_removed_edges:
                highlight_bonds_colors[bondIdx] = highlight_color_removed_dark
                highlight_bonds.append(bondIdx)
            elif (not modification_only) and (pair1 in result['removed_edges'] or pair2 in result['removed_edges']):
                highlight_bonds_colors[bondIdx] = highlight_color_removed
                highlight_bonds.append(bondIdx)

    img = draw_molecule(result['merged_mol'], highlightAtoms=highlight_atoms, highlightAtomColors=highlight_atoms_colors, highlightBonds=highlight_bonds, highlightBondColors=highlight_bonds_colors, output_type=output_type, **kwargs)
    if show_legend:
        legend = _generate_modification_legend(legend_font, output_type)
        img = _overlay_legend(img, legend, legend_position, output_type)
    return img


def draw_molecule_heatmap(mol, scores, output_type='png', show_labels = False, shrink_labels = False, annotation_scale = 1, show_legend = True, legend_width = 50, legend_font = 40, **kwargs):
    """
    Draw a molecule and color the atoms based on the scores
    
    Parameters
    ----------
    mol : rdkit molecule
        Molecule to Draw the heatmap for
    scores : list
        list of scores for each atom (the order should be the same as the order of atoms in the molecule)
    output_type : str, optional (default='png')
        type of output (png or svg)
    show_labels : bool, optional (default=False)
        add labels to the atoms
    shrink_labels : bool, optional (default=False)
        shrink the labels to more compact form
    annotation_scale : float, optional (default=1)
        size of the labels (font size)
    show_legend : bool, optional (default=True)
        show the legend or not
    legend_width : int, optional (default=50)
        width of the legend in pixels
    legend_font : int, optional (default=40)
        font size of the legend
        
    Returns
    -------
    img : numpy array or str
    
    Examples
    -------
    .. code-block:: python
    
        import modifinder.utilities.visualizer as mf_vis
        from matplotlib import pyplot as plt
        import numpy as np
        mol = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
        scores = np.random.rand(mol.GetNumAtoms())
        img = mf_vis.draw_molecule_heatmap(mol, scores, label="Caffeine", show_labels=True, legend_font=20)
        plt.imshow(img)
        plt.axis('off')
        plt.show()
        
    .. image:: ../_static/draw_molecule_heatmap.png
        :width: 300px
    """
    #TODO: test svg
    
    molecule = mu._get_molecule(mol)
    if type(molecule) == str or molecule is None:
        raise ValueError("Molecule not found")
    
    # get copy of the molecule
    mol = Chem.Mol(molecule)

    if min(scores) < 0:
        scores = [x + abs(min(scores)) for x in scores]
    if max(scores) > 0:
        vals = [x/max(scores) for x in scores]
    else:
        vals = [0 for x in scores]
    
    colors = dict()
    for i in range(0, mol.GetNumAtoms()):
        colors[i] = _get_heat_map_colors(vals[i])
    
    if show_labels:
        for atom in mol.GetAtoms():
            lbl = str(round(scores[atom.GetIdx()], 2))
            if shrink_labels:
                if scores[atom.GetIdx()] == 0:
                    lbl = ""
                else:
                    lbl = str(int(round(scores[atom.GetIdx()], 2)*100))
            atom.SetProp('atomNote',lbl)

    # # set the colors in kwargs
    kwargs['highlightAtoms'] = list(range(mol.GetNumAtoms()))
    kwargs['highlightAtomColors'] = colors
    kwargs['highlightBonds'] = []

    kwargs['extra_info'] = {}
    kwargs['extra_info']['annotation_scale'] = annotation_scale

    kwargs['extra_info']['output_type'] = output_type
    kwargs['extra_info']['show_legend'] = show_legend
    kwargs['extra_info']['legend_width'] = legend_width
    kwargs['extra_info']['legend_font'] = legend_font
    kwargs['extra_info']['scores'] = scores

    # draw the molecule
    img = draw_molecule(mol, output_type=output_type, **kwargs)
    return img


def draw_spectrum(spectrum, output_type='png', normalize_peaks = False, colors: dict = {}, flipped = False, size = None, show_x_label = False, show_y_label = False, font_size = None, bar_width = 3, x_lim = None,  dpi = 300, **kwargs):
    """
    Draw a spectrum 
    
    Parameters
    ----------
    spectrum : Spectrum or list of tuples or USI or Accession
        Spectrum to draw
    output_type : str, optional (default='png')
        type of output (png or svg)
    normalize_peaks : bool, optional (default=False)
        normalize the peaks or not
    colors : dict, optional (default={})
        dictionary of colors for the peaks, keys are the indices of the peaks,
        if value is a list, the first value is the color of top half and the 
        second value is the color of the bottom half.
    flipped : bool, optional (default=False)
        flip the spectrum or not
    size : tuple, optional (default=None)
        size of the figure
    show_x_label : bool, optional (default=False)
        show x label or not
    show_y_label : bool, optional (default=False)
        show y label or not
    font_size : int, optional (default=None)
        font size of the labels
    bar_width : int, optional (default=3)
        width of the spectrum bars
    x_lim : tuple, optional (default=None)
        x limit of the figure
    dpi : int, optional (default=300)
        dpi of the stored image
    kwargs : dict
        additional arguments for drawing the spectrum in matplotlib
    
    Returns
    -------
    img : numpy array or str
    
    Examples
    -------
    .. code-block:: python
    
        import modifinder.utilities.visualizer as mf_vis
        from matplotlib import pyplot as plt
        import numpy as np
        mz = [100, 200, 270, 400, 450]
        intensity = [0.1, 0.2, 0.5, 0.4, 0.3]
        colors = {
            0: 'red',
            1: ['blue', 'green'],
            3: '#FFA500',
            4: (0.9,0.9,0.2)
        }
        img = mf_vis.draw_spectrum(list(zip(mz, intensity)), output_type='png', show_x_label=True, show_y_label=True, colors=colors)
        plt.imshow(img)
        plt.axis('off')
        plt.show()
        
    .. image:: ../_static/draw_spectrum.png
        :width: 300px
    """
    spectrum = to_spectrum(spectrum)
    if normalize_peaks:
        spectrum.normalize_peaks()
    
    # if output type is ax, use the ax to draw the spectrum
    if isinstance(output_type, plt.Axes):
        ax = output_type
    else:
        if size is None:
            fig, ax = plt.subplots()
        else:
            fig, ax = plt.subplots(figsize=size)
    if normalize_peaks:
        ax.set_ylim(0, 1)
        ax.set_yticks(list(np.arange(0, 1.001, 0.2)))

    if font_size is not None:
        plt.rcParams.update({'font.size': font_size})

    for i in range(len(spectrum.mz)):
        mz = spectrum.mz[i]
        intensity = spectrum.intensity[i]
        color = colors.get(i, 'gray')
        if type(color) == list:
            if color[0] is None:
                color[0] = 'gray'
            if color[1] is None:
                color[1] = 'gray'
            
            ax.bar(mz, intensity/2, color=color[0], width=bar_width, bottom=intensity/2)
            ax.bar(mz, intensity/2, color=color[1], width=bar_width)
        else:
            ax.bar(mz, intensity, color=color, width=bar_width)
    
    if show_x_label:
        ax.set_xlabel("m/z")
    if show_y_label:
        ax.set_ylabel("Intensity")
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if flipped:
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.spines['top'].set_visible(True)
        ax.spines['bottom'].set_visible(False)
    
    if x_lim is not None:
        ax.set_xlim(x_lim)

    if output_type == "png":    
        alignment_bytes = BytesIO()
        fig.patch.set_alpha(0)  # Set the figure patch alpha to 0 (fully transparent)
        fig.savefig(alignment_bytes, format='png', bbox_inches='tight', pad_inches=0, transparent=True, dpi=dpi)
        alignment_bytes.seek(0)
        plt.close(fig)
        # convert the image to numpy array
        img = mpimg.imread(alignment_bytes)
        return img
    
    elif output_type == "svg":
        fig.patch.set_alpha(0)
        fig.canvas.draw()
        buffer = io.BytesIO()
        fig.savefig(buffer, format='svg', bbox_inches='tight', pad_inches=0, transparent=True)
        buffer.seek(0)
        svg = buffer.read().decode('utf-8')
        plt.close(fig)
        return svg
    elif output_type == "ax":
        return ax


def draw_alignment(spectrums, matches = None, output_type='png', normalize_peaks = False, size = None, dpi=300, draw_mapping_lines = True, ppm=40, x_lim=None, **kwargs):
    """
    Draw the alignment of multiple spectrums
    
    Parameters
    ----------
    spectrums : list of SpectrumTuple or list of list of tuples (mz, intensity)
        list of spectrums to draw
    matches : list of list of tuples or list of tuples, optional (default=None)
        matching between the spectrums
    output_type : str, optional (default='png')
        type of output (png or svg)
    normalize_peaks : bool, optional (default=False)
        normalize the peaks or not
    size : tuple, optional (default=None)
        size of the figure
    dpi : int, optional (default=300)
        dpi of the stored image
    draw_mapping_lines : bool, optional (default=True)
        draw the mapping lines
    ppm : float, optional (default=40)
        ppm value for the alignment
    x_lim : tuple, optional (default=None)
        x limit of the figure
    kwargs : dict
        additional arguments for drawing the spectrum in matplotlib
    
    Returns
    -------
    img : numpy array or str
    
    Examples
    -------
    .. code-block:: python
    
        import modifinder.utilities.visualizer as mf_vis
        from matplotlib import pyplot as plt
        import numpy as np
        peaks1 = list(zip([100, 200, 270, 400, 450], [0.1, 0.2, 0.5, 0.4, 0.3]))
        peaks2 = list(zip([100, 230, 350, 360, 430], [0.1, 0.3, 0.2, 0.4, 0.5]))
        peaks3 = list(zip([120, 230, 300, 380, 550], [0.1, 0.2, 0.5, 0.4, 0.3]))
        matches = [[(0, 0), (1, 1), (3, 4)], [(0, 0), (1, 1), (3, 3)]]
        img = mf_vis.draw_alignment([peaks1, peaks2, peaks3], matches=matches,
                                    output_type='png',
                                    normalize_peaks=True,
                                    x_lim=(0, 550))
        plt.imshow(img)
        plt.axis('off')
        plt.show()

    .. image:: ../_static/draw_alignment.png
        :width: 400px
    """
    
    spectrums = [to_spectrum(spectrum) for spectrum in spectrums]
    if normalize_peaks:
        for spectrum in spectrums:
            spectrum.normalize_peaks()

    if x_lim is None:
        x_lim = (min(spectrums[0].mz), max(spectrums[0].mz))
        for spectrum in spectrums:
            x_lim = (min(x_lim[0], min(spectrum.mz)), max(x_lim[1], max(spectrum.mz)))
    
    # calculating the colors
    # matches should be a list of matchings for each spectrum pairs, each matching is a list of tuples (index1, index2)

    colors = [dict() for i in range(len(spectrums))]

    if matches is None:
        # perform the alignment
        matches = []
    
    # if matches == 'default':
    #     matches = []
    #     for i in range(len(spectrums) - 1):
    #         cosine, match = _cosine_fast(spectrums[i], spectrums[i+1], 0.1, ppm, True)
    #         matches.append(match)
            
    
    if len(matches) > 0 and len(matches) != len(spectrums) - 1:
        raise ValueError("Number of matches should be equal to the number of spectrums - 1")
    
    # in case there is only one match, accept it as both a list of tuples or a list of lists
    if len(matches) == 1:
        if type(matches[0][0]) == int:
            matches = [matches]

    flipped = False
    if "flipped" in kwargs:
        flipped = kwargs["flipped"]
        kwargs.pop("flipped")

    lines = []
    for match_index, match in enumerate(matches):
        for pair in match:
            temp_color = 'blue'
            if is_shifted(spectrums[match_index].mz[pair[0]], spectrums[match_index+1].mz[pair[1]], ppm, None):
                temp_color = 'red'
            
            if pair[0] not in colors[match_index]:
                if match_index > 0:
                    colors[match_index][pair[0]] = [None, temp_color]
                else:
                    colors[match_index][pair[0]] = temp_color
            else:
                colors[match_index][pair[0]][1] = temp_color
            
            if pair[1] not in colors[match_index+1]:
                if match_index < len(spectrums) - 2:
                    colors[match_index+1][pair[1]] = [temp_color, None]
                else:
                    colors[match_index+1][pair[1]] = temp_color
            else:
                colors[match_index+1][pair[1]][0] = temp_color
            
            if draw_mapping_lines:
                if flipped and match_index == 0:
                    lines.append([(match_index, spectrums[match_index].mz[pair[0]], 0), (match_index+1, spectrums[match_index+1].mz[pair[1]], 0), temp_color])
                else:
                    lines.append([(match_index, spectrums[match_index].mz[pair[0]], 0), (match_index+1, spectrums[match_index+1].mz[pair[1]], spectrums[match_index+1].intensity[pair[1]]), temp_color])
    
                
    if size is None:
        fig, axs = plt.subplots(len(spectrums), 1)
    else:
        fig, axs = plt.subplots(len(spectrums), 1, figsize=size)

    for index, spectrum in enumerate(spectrums):
        if len(spectrums) == 2 and index == 1:
            draw_spectrum(spectrum, output_type=axs[index], colors=colors[index], x_lim=x_lim, normalize_peaks=normalize_peaks, flipped=flipped, **kwargs)
        else:
            draw_spectrum(spectrum, output_type=axs[index], colors=colors[index], x_lim=x_lim, normalize_peaks=normalize_peaks, flipped=False, **kwargs)

    
    for line in lines:
        if "bar_width" in kwargs:
            bar_width = max(1, kwargs["bar_width"]//1.5)
            con = ConnectionPatch(xyA=(line[1][1], line[1][2]), xyB=(line[0][1], line[0][2]), coordsA="data", coordsB="data", axesA=axs[line[1][0]], axesB=axs[line[0][0]], color=line[2], linestyle='dotted', linewidth=bar_width)
        else:
            con = ConnectionPatch(xyA=(line[1][1], line[1][2]), xyB=(line[0][1], line[0][2]), coordsA="data", coordsB="data", axesA=axs[line[1][0]], axesB=axs[line[0][0]], color=line[2], linestyle='dotted')
        axs[line[1][0]].add_artist(con)
    
    # remove all extra paddings (but make sure the labels are not cut)
    if flipped and len(spectrums) == 2:
        # remove xtick labels
        axs[1].set_xticklabels([])
        plt.subplots_adjust(hspace=0.01, wspace=0)
    else:
        plt.subplots_adjust(hspace=0.1, wspace=0)
    plt.tight_layout()
    
    if output_type == "png":
        alignment_bytes = BytesIO()
        fig.patch.set_alpha(0)  # Set the figure patch alpha to 0 (fully transparent)
        fig.savefig(alignment_bytes, format='png', bbox_inches='tight', pad_inches=0, transparent=True, dpi=dpi)
        alignment_bytes.seek(0)
        plt.close(fig)
        # convert the image to numpy array
        img = mpimg.imread(alignment_bytes)
        return img
    elif output_type == "svg":
        fig.patch.set_alpha(0)
        fig.canvas.draw()
        buffer = io.BytesIO()
        fig.savefig(buffer, format='svg', bbox_inches='tight', pad_inches=0, transparent=True)
        buffer.seek(0)
        svg = buffer.read().decode('utf-8')
        plt.close(fig)
        return svg
    else:
        return fig, axs


def draw_frag_of_molecule(mol, fragment, output_type='png', **kwargs):
    """
    draws a fragment of the molecule. The fragment is represented by a binary string where 1 indicates the presence of the atom and 0 indicates the absence.
    
    Parameters
    ----------
    mol : rdkit molecule
        Molecule to draw the fragment for
    fragment : int or list
        fragment represented in base 10 where in the binary representation 1 indicates the presence of the atom and 0 indicates the absence, or a list of atom indices
    output_type : str, optional (default='png')
        type of output (png or svg)
    kwargs : dict
        additional arguments for drawing the molecule in rdkit
        
    Returns
    -------
    img : numpy array or str
    
    Examples
    -------
    .. code-block:: python
    
        import modifinder.utilities.visualizer as mf_vis
        from matplotlib import pyplot as plt
        from rdkit import Chem
        mol = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
        def mol_with_atom_index(mol):
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx())
            return mol
        mol = mol_with_atom_index(mol)
        fragment = int("110111", 2) # Convert binary to decimal
        img = mf_vis.draw_frag_of_molecule(mol, fragment, output_type='png')
        plt.imshow(img)
        plt.axis('off')
        plt.show()
        
    .. image:: ../_static/draw_frag_of_molecule.png
        :width: 300px
        
    """
    print("in visualizer", fragment)
    mol = mu._get_molecule(mol)
    highlightAtoms = []
    if isinstance(fragment, int):
        for i in range(0, mol.GetNumAtoms()):
            if fragment & (1 << i):
                highlightAtoms.append(i)
    elif isinstance(fragment, list):
        highlightAtoms = fragment
    # print(highlightAtoms)
    kwargs['highlightAtoms'] = highlightAtoms
    img = draw_molecule(mol, output_type=output_type, **kwargs)
    print(img)
    return img


def _get_heat_map_colors(val):
    """
    Get the heat map color based on an input value
    :param val: float - value between 0 and 1
    """
    if val == 0:
        return (val, 0.2*(1-val), (1-val), 0.4)
    else:
        return (0.95*val, 0.2*(1-val), (1-val), val*0.3 + 0.40)


def _add_heatmap_legend(img, scores, colors, image_width, image_height, legend_span, fontSize, output_type):
    """
    Add legend to the image
    :param img: image
    :param scores: list of scores
    :param colors: dictionary of colors
    :param type: str - type of output (png or svg)
    """
    steps = 100
    if output_type == "png":
        legend_patch = np.zeros((legend_span, image_width, 4))
        for i in range(0, steps):
                legend_patch[:, int(i*image_width/steps):int((i+1)*image_width/steps), :] = np.array([_get_heat_map_colors(i/steps)]*legend_span).reshape(legend_span, 1, 4)
    
        legend_patch[:,:,:] = legend_patch[:,:,:]*255
        legend_image = Image.fromarray(legend_patch.astype(np.uint8), 'RGBA')

        # add legend text
        draw = ImageDraw.Draw(legend_image)
        path_of_this_file = os.path.abspath(os.path.dirname(__file__))
        font = ImageFont.truetype(os.path.join(path_of_this_file, "fonts/NotoSans-Regular.ttf"), fontSize)
        draw.text((5, legend_span//2), "Low likelihood", (255, 255, 255), font=font, anchor="lm")
        draw.text((image_width-5, legend_span//2), "high likelihood", (255, 255, 255), font=font, anchor="rm")

        legend_patch = np.array(legend_image)

        # if image is in RGB format, convert it to RGBA
        if img.max() <= 1:
            img = (img*255).astype(np.uint8)
        if img.shape[2] == 3:
            img = np.concatenate((img, np.ones((img.shape[0], img.shape[1], 1), dtype=np.uint8)*255), axis=2)
        img = np.concatenate((img, legend_patch), axis=0)
        return img
    else:
        svg = """<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg">""".format(legend_span, image_width)
        rectWidth = legend_span/steps
        rectHeight = image_width

        for i in range(steps, -1, -1):
            color = _get_heat_map_colors(i/steps)
            svg += """<rect x="{}" y="{}" width="{}" height="{}" fill="rgba({}, {}, {}, {})"/>""".format(i*rectWidth, image_height-legend_span, rectWidth*1.01, rectHeight,
                                                                                                          color[0]*255, color[1]*255,color[2]*255, color[3])
            svg += """<text x="{}" y="{}" font-size="{}px" fill="white">{}</text>""".format(5, image_height-legend_span + (rectHeight+fontSize/2)/2, fontSize, "Low likelihood")
            text = "high likelihood"
            # add text to the right side
            text = text[::-1]
            svg += """<text x="{}" y="{}" font-size="{}px" fill="white", text-anchor="end">{}</text>""".format(legend_span-5, image_height-legend_span + (rectHeight+fontSize/2)/2, fontSize, text)
        svg += "</svg>"
        return svg
    

def _generate_modification_legend(fontSize, output_type, **kwargs):
    """
    generate legend for image
    :param img: image
    :param type: str - type of output (png or svg)
    return: numpy array of the image

    """
    if output_type == "png":
        fig, ax = plt.subplots(figsize=(1, 1))
        categories = ['Common', 'Added', 'Removed']
        colors = [(0.1, 0.1, 0.8, 0.7), (0.1, 0.8, 0.1, 0.7), (0.8, 0.1, 0.1, 0.7)]

        # Plot the dummy data to create a legend
        for color, category in zip(colors, categories):
            ax.plot([], [], color=color, label=category, linewidth=fontSize//1.5)
        legend_location = (1, 1)
        legend = ax.legend(fontsize=fontSize, **kwargs, loc='center', bbox_to_anchor=legend_location)
        ax.axis('off')
        fig.patch.set_alpha(0)
        fig.canvas.draw()
        legend_bytes = BytesIO()
        fig.savefig(legend_bytes, format='png', bbox_inches='tight', pad_inches=0.1, transparent=True)
        legend_bytes.seek(0)
        plt.close(fig)
        legend_image = Image.open(legend_bytes)
        legend_image = np.array(legend_image)
        return legend_image
    else:
        # return empty svg
        return "<svg></svg>"
        # TODO: Implement this
        raise NotImplementedError("SVG not implemented")


def _overlay_legend(img, legend, position, output_type):
    """
    Overlay the legend on the image
    :param img: image
    :param legend: image
    :param position: tuple - position of the legend
    :param output_type: str - type of output (png or svg)
    """
    if output_type == "png":
        if img.max() <= 1:
            img = (img*255).astype(np.uint8)
        else:
            img = img.astype(np.uint8)
        img2 = Image.fromarray(img)
        legend = Image.fromarray(legend)

        # if position is not provided, place the legend at the bottom right corner
        if position is None:
            position = (max(0, img2.size[0] - legend.size[0]), max(img2.size[1] - legend.size[1], 0))
        
        if type(position) != tuple:
            raise ValueError("Position should be a tuple")
        
        size = (max(img2.size[0], legend.size[0] + position[0]), max(img2.size[1], legend.size[1] + position[1] ))
        combined_image = Image.new("RGBA", size, img2.getpixel((0, 0)))
        combined_image.paste(img2, (0, 0))
        combined_image.paste(legend, position, legend)
        return np.array(combined_image)
    else:
        return img
        # TODO: Implement this
        raise NotImplementedError("SVG not implemented")

def _cname2hex(cname):
    colors = dict(mpl.colors.BASE_COLORS, **mpl.colors.CSS4_COLORS) # dictionary. key: names, values: hex codes
    try:
       hex = colors[cname]
       return hex
    except KeyError:
       print(cname, ' is not registered as default colors by matplotlib!')
    return None

def _hex2rgb(hex, normalize=False):
    h = hex.strip('#')
    rgb = np.asarray(list(int(h[i:i + 2], 16) for i in (0, 2, 4)))
    return rgb

def _handle_color(color):
    if type(color) == str:
        if color[0] == '#':
            return _hex2rgb(color)
        else:
            return _hex2rgb(_cname2hex(color))
    else:
        return color
    

def return_public_functions():
    """
    Return the public functions of the module
    """
    return {
        "draw_molecule": draw_molecule,
        # "draw_modifications": draw_modifications,
        "draw_molecule_heatmap": draw_molecule_heatmap,
        # "draw_spectrum": draw_spectrum,
        # "draw_alignment": draw_alignment,
        "draw_frag_of_molecule": draw_frag_of_molecule
    }