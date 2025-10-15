import csv
import svg_functions
from typing import IO


def create_horizontal_axis(svg_output: IO, x1: int, y1: int, line_width: int, y2: int):
    """Draw the horizontal axis of coordinates.

    :param svg_output: svg file to write to
    :param x1: x coordinate of the axis start
    :param y1: y coordinate of the axis start
    :param line_width: length of the line
    :param y2: y coordinate of the axis end
    """
    # horizontal line
    svg_output.write(
        '<line x1="' + str(x1) + '" y1="' + str(y1) + '" x2="' + str(x1 + line_width) + '" y2="' + str(y2) +
        '" style="stroke: ' + 'black' + '; stroke-width:' + '1' + '; stroke-linecap: round" />')


def create_axis_labels(svg_output: IO, x1: int, y1: int, tick_height: int, pos: int, axis_start: int, axis_label: str):
    """Draw axis labels on the horizontal axis.

    :param svg_output: svg file to write to
    :param x1: x coordinate of the axis label
    :param y1: y coordinate of the axis label
    :param tick_height: axis tick height
    :param pos: mtDNA position for the label
    :param axis_start: mtDNA position corresponding to the start of the horizontal axis
    :param axis_label: the text label
    """
    # vertical lines and labels
    svg_output.write(  # axis tick
        '<line x1="' + str((pos - axis_start) + x1) + '" y1="' + str(y1) +
        '" x2="' + str((pos - axis_start) + x1) + '" y2="' + str(y1 + tick_height) +
        '" style="stroke: ' + 'black' + '; stroke-width:' + '1' + '; stroke-linecap: round" />')
    svg_output.write(  # axis label
        '<text text-anchor="middle" alignment-baseline="middle" '
        'x="' + str((pos - axis_start) + x1) + '" y="' + str(y1 + tick_height + 7.5) +
        '" style = "font-size:' + str(9.5) +
        '; fill: #000000; font-family: arial; font-weight: normal;" > ' + str('m.' + str(axis_label)) + '</text>')
    

def control_region(svg_output):
    """Visualize constraint for each non-coding element in the control region as a svg.

    """
    # print axis line and labels
    # horizontal line
    create_horizontal_axis(svg_output=svg_output, x1=50, y1=20, line_width=((600 + 16569) - 16000), y2=20)
    
    # vertical lines and labels
    for pos in [16000, 16100, 16200, 16300, 16400, 16500, 100, 200, 300, 400, 500, 600]:
        new_pos = (pos + 16569) if pos <= 600 else pos  # to handle artificial break
        create_axis_labels(
            svg_output=svg_output, x1=50, y1=20, tick_height=7, pos=new_pos, axis_start=16000, axis_label=str(pos))
    
    # manually draw axis label at m.16569/m.1
    svg_output.write(
        '<line x1="' + str((16569 - 16000) + 50) + '" y1="' + str(20) +
        '" x2="' + str((16569 - 16000) + 50) + '" y2="' + str(27.5) +
        '" style="stroke: ' + 'black' + '; stroke-width:' + '1.25' + '; stroke-linecap: round" />')
    svg_output.write(
        '<text text-anchor="middle" alignment-baseline="middle"  x="' + str((16569 - 16000) + 50) + '" y="' + str(35) +
        '" style = "font-size:' + str(9.5) +
        '; fill: #000000; font-family: arial; font-weight: normal;" > ' + str('m.16569/m.1') + '</text>')
        
    for row in csv.DictReader(open('output_files/oe/noncoding_obs_exp.txt'), delimiter='\t'):
        row1 = ['MT-HV1', 'MT-HV2', 'MT-HV3']
        row2 = ['MT-TAS2', 'MT-TAS', 'MT-CSB1', 'MT-CSB2', 'MT-CSB3', 'MT-LSP', 'MT-HSP1']
        row3 = ['MT-5', 'MT-3L', 'MT-TFX', 'MT-TFY', 'MT-HPR', 'MT-4H', 'MT-3H', 'MT-TFL', 'MT-TFH']
        row4 = ['']  # ['MT-OHR']
        row5 = ['']  # ['MT-7SDNA']
        # set coordinates
        x = ((int(row["start"]) - 16000) if (int(row["start"]) > 16000) else ((int(row["start"]) + 16569) - 16000)) + 50
        if row["locus"] in row1:
            y = 60
        elif row["locus"] in row2:
            y = 90
        elif row["locus"] in row3:
            y = 120
        elif row["locus"] in row4:
            y = 150
        elif row["locus"] in row5:
            y = 180
        else:
            y = ''
        width = (int(row["end"]) - int(row["start"])) if (int(row["start"]) < int(row["end"])) else \
            ((int(row["end"]) + 16569) - int(row["start"]))
        font_size = 9.5 if row["locus"] in (row1 + row2 + row3) else 9.5
        
        # fill color by OEUF
        fill_color = svg_functions.return_color(value=float(row["upper_CI"]))
        
        # set fill color and label for each locus
        if (row["locus"] in (row1 + row2 + row3 + row4 + row5)) and (width > 10):
            svg_output.write(
                '<rect x="' + str(x) + '" y="' + str(y) + '" width="' + str(width) + '" height="20" style="fill:' +
                fill_color + ';stroke:rgb(0,0,0);stroke-width:1" />')
            
            legend = '<text text-anchor="start" x="' + str(x + width) + '" y="' + str(y + 12.5) + \
                     '" style = "font-size:' + str(font_size) + \
                     '; fill: #000000; font-family: arial; font-weight: normal;" > ' + \
                     str('-' + row["locus"].replace('MT-', '')) + '</text>'
            
            # these loci need the label on the left
            if row["locus"] in ['MT-TAS2', 'MT-TFX', 'MT-4H']:
                legend = '<text text-anchor="end" x="' + str(x) + '" y="' + str(y + 12.5) + \
                         '" style = "font-size:' + str(font_size) + \
                         '; fill: #000000; font-family: arial; font-weight: normal;" > ' + \
                         str(row["locus"].replace('MT-', '') + '-') + '</text>'
                
            svg_output.write(legend)

    svg_functions.color_legend(
        svg_output=svg_output, x=((600 + 16569) - 16000 + 100), y1=48, width=25, height=90, font_size=10)


def control_region_plain(svg_output):
    """Visualize the position of the non-coding elements in the control region.

    """
    # print axis line and labels
    # horizontal line
    create_horizontal_axis(svg_output=svg_output, x1=50, y1=20, line_width=((600 + 16569) - 16000), y2=20)
    
    # vertical lines and labels
    for pos in [16000, 16100, 16200, 16300, 16400, 16500, 100, 200, 300, 400, 500, 600]:
        new_pos = (pos + 16569) if pos <= 600 else pos  # to handle artificial break
        create_axis_labels(
            svg_output=svg_output, x1=50, y1=20, tick_height=7, pos=new_pos, axis_start=16000, axis_label=str(pos))
    
    # manually draw axis label at m.16569/m.1
    svg_output.write(
        '<line x1="' + str((16569 - 16000) + 50) + '" y1="' + str(20) +
        '" x2="' + str((16569 - 16000) + 50) + '" y2="' + str(27.5) +
        '" style="stroke: ' + 'black' + '; stroke-width:' + '1.25' + '; stroke-linecap: round" />')
    
    for row in csv.DictReader(open('output_files/oe/noncoding_obs_exp.txt'), delimiter='\t'):
        row1 = ['MT-HV1', 'MT-HV2', 'MT-HV3']
        row2 = ['MT-TAS2', 'MT-TAS', 'MT-CSB1', 'MT-CSB2', 'MT-CSB3', 'MT-LSP', 'MT-HSP1']
        row3 = ['MT-5', 'MT-3L', 'MT-TFX', 'MT-TFY', 'MT-HPR', 'MT-4H', 'MT-3H', 'MT-TFL', 'MT-TFH']
        row4 = ['']  # ['MT-OHR']
        row5 = ['']  # ['MT-7SDNA']
        # set coordinates
        x = ((int(row["start"]) - 16000) if (int(row["start"]) > 16000) else ((int(row["start"]) + 16569) - 16000)) + 50
        if row["locus"] in row1:
            y = 60
        elif row["locus"] in row2:
            y = 90
        elif row["locus"] in row3:
            y = 120
        elif row["locus"] in row4:
            y = 150
        elif row["locus"] in row5:
            y = 180
        else:
            y = ''
        width = (int(row["end"]) - int(row["start"])) if (int(row["start"]) < int(row["end"])) else \
            ((int(row["end"]) + 16569) - int(row["start"]))
        font_size = 9.5 if row["locus"] in (row1 + row2 + row3) else 9.5
        
        # fill color - set to be grey
        fill_color = "#BEBEBE"

        # set fill color and label for each locus
        if (row["locus"] in (row1 + row2 + row3 + row4 + row5)) and (width > 10):
            svg_output.write(
                '<rect x="' + str(x) + '" y="' + str(y) + '" width="' + str(width) + '" height="20" style="fill:' +
                fill_color + ';stroke:rgb(0,0,0);stroke-width:1" />')
            
            legend = '<text text-anchor="start" x="' + str(x + width) + '" y="' + str(y + 12.5) + \
                     '" style = "font-size:' + str(font_size) + \
                     '; fill: #000000; font-family: arial; font-weight: normal;" > ' + \
                     str('-' + row["locus"].replace('MT-', '')) + '</text>'
            
            # these loci need the label on the left
            if row["locus"] in ['MT-TAS2', 'MT-TFX', 'MT-4H']:
                legend = '<text text-anchor="end" x="' + str(x) + '" y="' + str(y + 12.5) + \
                         '" style = "font-size:' + str(font_size) + \
                         '; fill: #000000; font-family: arial; font-weight: normal;" > ' + \
                         str(row["locus"].replace('MT-', '') + '-') + '</text>'
                
            svg_output.write(legend)


def draw_figures():
    """Visualize non-coding control region constraint as a svg.
    
    """
    with open('figure_scripts/figures/FigureED7_control_region.svg', 'w') as svg_output:
        svg_output.write('<svg height="200" width="1400" xmlns="http://www.w3.org/2000/svg">\n\n')
        svg_output.write('\n<!-- MARKERS -->\n')
        
        control_region(svg_output)
        
        svg_output.write('</svg>')

    with open('figure_scripts/extended_data_figures/FigureED8a_control_region.svg', 'w') as svg_output:
        svg_output.write('<svg height="200" width="1300" xmlns="http://www.w3.org/2000/svg">\n\n')
        svg_output.write('\n<!-- MARKERS -->\n')
    
        control_region_plain(svg_output)
    
        svg_output.write('</svg>')


if __name__ == "__main__":
    draw_figures()
