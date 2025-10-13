import csv
import os


def pair_info():
    """Generate dictionary that will be used to draw rRNA svgs.

    :return: matches, dictionary with information on bases in pairs, and base coordinates
    """
    matches = {}
    for file in os.listdir('required_files/svg_input/'):
        if file.endswith(".tsv") and ('RNR' in file):
            gene = file.split('_')[0]
            
            for row in csv.DictReader(open('required_files/svg_input' + '/' + file), delimiter='\t'):
                infoType = row["Type"]
                genCoord = row["Genomic_coordinate"]
                base = row["RNA.base"]
                x1 = row["x1"]
                y1 = row["y1"]
                pair = row["Pair_coordinate"]
                
                if infoType == "b" and pair:
                    matches[(pair, gene)] = (genCoord, base, x1, y1)
    return matches


def draw_rrna_svgs(matches):
    """Visualize rRNA regional constraint as a svg.

    :param matches: dictionary with information on bases in pairs, and base coordinates
    """
    rc_dict = {}
    for row in csv.DictReader(
            open('output_files/regional_constraint/mito_regional_constraint_annotation.txt'), delimiter='\t'):
        rc_dict[row["POS"]] = row["in_rc"]

    for gene in ['MT-RNR1', 'MT-RNR2']:
        file_name = 'extended_data_figures/FigureED5g_RNR1' if gene == "MT-RNR1" else 'extended_data_figures/FigureED5h_RNR2'
        with open('figure_scripts/%s.svg' % file_name, 'w') as svg_output:
            svg_output.write('<svg height="2000" width="3000" xmlns="http://www.w3.org/2000/svg">\n\n')
            svg_output.write('\n<!-- MARKERS -->\n')

            for row in csv.DictReader(
                    open('required_files/svg_input/%s_with_pairs.tsv' % gene), delimiter='\t'):
                infoType = row["Type"]
                genCoord = row["Genomic_coordinate"]
                base = row["RNA.base"]
                x1 = row["x1"]
                y1 = row["y1"]
                pair = row["Pair_coordinate"]
                font = "monospace"
                fontWeight = 'normal'
                fontSize = '15'
                title = genCoord
                lineColor = "#000000"
                strokeWidth = '1'
                circleColor = "#000000"
                radius = "2"  # circle radius

                if gene == 'MT-RNR1':
                    fontSize = '8'
                    strokeWidth = '0.7'
                    radius = "1.3"
                    # below are to draw RNR1 non-WC lines
                    svg_output.write(
                        '<line x1="' + str(486) + '" y1="' + str(410) + '" x2="' + str(494) + '" y2="' + str(420) + '" style="stroke: ' + lineColor + '; stroke-width:1 ; stroke-linecap: round" >' + '<title>' + str(1137) + ',' + str(1138) + '</title> </line>')
                    svg_output.write(
                        '<line x1="' + str(540) + '" y1="' + str(503) + '" x2="' + str(529) + '" y2="' + str(509) + '" style="stroke: ' + lineColor + '; stroke-width:1 ; stroke-linecap: round" >' + '<title>' + str(657) + ',' + str(658) + '</title> </line>')
                    svg_output.write(
                        '<line x1="' + str(504) + '" y1="' + str(508) + '" x2="' + str(520) + '" y2="' + str(511) + '" style="stroke: ' + lineColor + '; stroke-width:1 ; stroke-linecap: round" >' + '<title>' + str(656) + ',' + str(657) + '</title> </line>')
                    svg_output.write(
                        '<line x1="' + str(316) + '" y1="' + str(606) + '" x2="' + str(315) + '" y2="' + str(623) + '" style="stroke: ' + lineColor + '; stroke-width:1 ; stroke-linecap: round" >' + '<title>' + str(687) + ',' + str(688) + '</title> </line>')
                    svg_output.write(
                        '<line x1="' + str(315) + '" y1="' + str(635) + '" x2="' + str(315) + '" y2="' + str(648) + '" style="stroke: ' + lineColor + '; stroke-width:1 ; stroke-linecap: round" >' + '<title>' + str(688) + ',' + str(689) + '</title> </line>')
                    svg_output.write(
                        '<line x1="' + str(235) + '" y1="' + str(775) + '" x2="' + str(238) + '" y2="' + str(755) + '" style="stroke: ' + lineColor + '; stroke-width:1 ; stroke-linecap: round" >' + '<title>' + str(802) + ',' + str(801) + '</title> </line>')
                    svg_output.write(
                        '<line x1="' + str(236) + '" y1="' + str(785) + '" x2="' + str(242) + '" y2="' + str(802) + '" style="stroke: ' + lineColor + '; stroke-width:1 ; stroke-linecap: round" >' + '<title>' + str(801) + ',' + str(800) + '</title> </line>')
                    svg_output.write(
                        '<path d="M426,462 Q415,443 426,426" stroke="black" fill="transparent" />')

                if gene == 'MT-RNR2':   # draw the non-WC line
                    svg_output.write(
                        '<line x1="' + str(1170) + '" y1="' + str(1018) + '" x2="' + str(1762) + '" y2="' + str(1018) + '" style="stroke: ' + lineColor + '; stroke-width:2 ; stroke-linecap: round" >' + '<title>' + str(2452) + ',' + str(2453) + '</title> </line>')
                
                # set coloring of regionally constrained, modified and disease-associated bases
                if genCoord:
                    if int(genCoord) in [1076, 1486, 1488, 1583, 1584, 2617, 2815, 3039, 3040, 3067]:  # modified bases
                        textColor = "blue"
                        fontWeight = 'bold'
                        fontSize = '16' if gene == "MT-RNR1" else '23'
                    elif int(genCoord) in [1494, 1555]:  # disease variant bases
                        textColor = "#7f4188"
                        fontWeight = 'bold'
                        fontSize = '16' if gene == "MT-RNR1" else '23'
                    elif rc_dict[genCoord] == "yes":
                        textColor = "#ff4040"
                        fontWeight = 'bold'
                        fontSize = '13' if gene == "MT-RNR1" else '20'
                    else:
                        textColor = "#984ea3"
                else:
                    if not genCoord:
                        textColor = "dimgrey"
                    else:
                        textColor = "#000000"

                base_svg = '<text  text-anchor="middle" alignment-baseline="middle" x="' + x1 + '" y="' + y1 + '" style = "font-size: ' + fontSize + '; fill: ' + textColor + '; font-family: ' + font + '; font-weight: ' + fontWeight + ';" >' + base + '<title>' + title + '</title> </text>'

                if infoType == "b":  # input only contains 'b' - except for 'n' at 3107 which is N
                    if pair:
                        if (genCoord, gene) in matches:
                            svg_output.write(base_svg)
                            base1 = base
                            genCoord1 = genCoord
                            x1 = int(x1)
                            y1 = int(y1)

                            genCoord2 = matches[(genCoord, gene)][0]
                            base2 = matches[(genCoord, gene)][1]
                            x2 = int(matches[(genCoord, gene)][2])
                            y2 = int(matches[(genCoord, gene)][3])

                            if x1 == x2:
                                middle_x = float((x1 + x2) / 2)
                                x3 = x1
                                x4 = x2

                                middle_y = float((y1 + y2) / 2)
                                y3 = str((middle_y + 4))
                                y4 = str((middle_y - 4))
                                
                            elif y1 == y2:
                                middle_x = float((x1 + x2) / 2)
                                x3 = str((middle_x + 4))
                                x4 = str((middle_x - 4))

                                middle_y = float((y1 + y2) / 2)
                                y3 = y1
                                y4 = y2

                            else:
                                
                                middle_x = float((x1 + x2) / 2)
                                length_x = abs(float((x1 - x2)) / 3)

                                middle_y = float((y1 + y2) / 2)
                                length_y = abs(float((y1 - y2)) / 3)

                                if (x1 > x2) and (y1 < y2):

                                    x3 = str(float((x1 - length_x)))
                                    y3 = str(float((y1 + length_y)))

                                    x4 = str(float((x2 + length_x)))
                                    y4 = str(float((y2 - length_y)))

                                elif (x1 < x2) and (y1 < y2):
                                    
                                    x3 = str(float((x1 + length_x)))
                                    y3 = str(float((y1 + length_y)))

                                    x4 = str(float((x2 - length_x)))
                                    y4 = str(float((y2 - length_y)))

                                elif (x1 > x2) and (y1 > y2):

                                    x3 = str(float((x1 - length_x)))
                                    y3 = str(float((y1 - length_y)))

                                    x4 = str(float((x2 + length_x)))
                                    y4 = str(float((y2 + length_y)))

                                elif (x1 < x2) and (y1 > y2):
                                    
                                    x3 = str(float((x1 + length_x)))
                                    y3 = str(float((y1 - length_y)))

                                    x4 = str(float((x2 - length_x)))
                                    y4 = str(float((y2 + length_y)))

                            svg_output.write(base_svg)
                            line = '<line x1="' + str(x3) + '" y1="' + str(y3) + '" x2="' + str(x4) + '" y2="' + str(y4) + '" style="stroke: ' + lineColor + '; stroke-width:' + strokeWidth + '; stroke-linecap: round" >' + '<title>' + genCoord1 + ',' + genCoord2 + '</title> </line>'
                            circle = '<circle cx="' + str(middle_x) + '" cy="' + str(middle_y) + '" r="' + radius + '" style="fill: ' + circleColor + ';" >' + '<title>' + genCoord1 + ',' + genCoord2 + '</title> </circle>'

                            if base1 == 'A':
                                if base2 == 'T':
                                    svg_output.write(line)
                                else:
                                    svg_output.write(circle)
                            if base1 == 'T':
                                if base2 == 'A':
                                    svg_output.write(line)
                                else:
                                    svg_output.write(circle)
                            if base1 == 'C':
                                if base2 == 'G':
                                    svg_output.write(line)
                                else:
                                    svg_output.write(circle)
                            if base1 == 'G':
                                if base2 == 'C':
                                    svg_output.write(line)
                                else:
                                    svg_output.write(circle)
                    else:
                        svg_output.write(base_svg)
                        
            if gene == "MT-RNR1":
                svg_output.write(
                    '<rect x="540" y="535" width="180" height="120" style = "fill:none;stroke-width:1.5;stroke:rgb(0,0,0)"> </rect>')
                
            svg_output.write('</svg>')


if __name__ == "__main__":
    matches = pair_info()
    draw_rrna_svgs(matches)
