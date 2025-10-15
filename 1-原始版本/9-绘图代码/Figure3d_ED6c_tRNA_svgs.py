import csv
import svg_functions


def draw_trna_pos_svg():
    """Visualize tRNA positional constraint as a svg.
    
    """
    # using MT-TL1 as framework, so using genomic and pair coordinates here
    # first get information to draw lines between pairs
    matches_pairs = {}
    for row in csv.DictReader(open('required_files/svg_input/tRNA_model.tsv'), delimiter='\t'):
        if row["Pair_coordinate"]:
            matches_pairs[row["Pair_coordinate"]] = (row["x1"], row["y1"])
    
    # now build dictionary of site constraint, using OEUF
    matches_site = {}
    for row in csv.DictReader(open('output_files/oe/tRNA_position_obs_exp.txt'), delimiter='\t'):
        matches_site[row["tRNA_position"]] = row["upper_CI"]

    # now draw the svg
    with open('figure_scripts/figures/Figure3d_tRNA_pos.svg', 'w') as svg_output:
        svg_output.write('<svg height="350" width="350" xmlns="http://www.w3.org/2000/svg">\n\n')
        svg_output.write('\n<!-- MARKERS -->\n')

        for row in csv.DictReader(open('required_files/svg_input/tRNA_model.tsv'), delimiter='\t'):
            site = row["Site_number"]
            x1 = row["x1"]
            y1 = row["y1"]
            font = "arial"
            fontWeight = 'normal'
            fontSize = '7'
            lineColor = "#000000"
            strokeWidth = '1'
            radius = "5.5"  # circle radius
            textColor = "#000000"

            circleColor = svg_functions.return_color(value=float(matches_site[site]))

            if ":a" in site:
                fontSize = '6'
                site = site.replace(':', "")
                
            circle = '<circle cx="' + str(x1) + '" cy="' + str(y1) + '" r="' + \
                     radius + '" style="fill: ' + circleColor + ';stroke:rgb(0,0,0);stroke-width:0.2" />'
            svg_output.write(circle)
            
            base_svg = '<text  text-anchor="middle" alignment-baseline="middle" x="' + x1 + '" y="' + y1 + \
                       '" style = "font-size: ' + fontSize + '; fill: ' + textColor + '; font-family: ' + font + \
                       '; font-weight: ' + fontWeight + ';" >' + site + '</text>'

            if row["Pair_coordinate"]:
                svg_output.write(base_svg)
                x1 = int(x1)
                x2 = int(matches_pairs[row["Genomic_coordinate"]][0])
                y1 = int(y1)
                y2 = int(matches_pairs[row["Genomic_coordinate"]][1])

                if x1 == x2:
                    x3 = x1
                    x4 = x2
                    middle_y = float((y1 + y2) / 2)
                    y3 = str((middle_y + 4))
                    y4 = str((middle_y - 4))
                elif y1 == y2:
                    middle_x = float((x1 + x2) / 2)
                    x3 = str((middle_x + 4))
                    x4 = str((middle_x - 4))
                    y3 = y1
                    y4 = y2
                else:
                    length_x = abs(float((x1 - x2)) / 3)
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
                        
                line = '<line x1="' + str(x3) + '" y1="' + str(y3) + '" x2="' + str(x4) + '" y2="' + str(y4) + \
                       '" style="stroke: ' + lineColor + '; stroke-width:' + strokeWidth + '; stroke-linecap: round" />'
                svg_output.write(line)  # all as lines

            else:
                svg_output.write(base_svg)

        svg_functions.color_legend(svg_output=svg_output, x=55, y1=80, width=20, height=65, font_size=7)

        svg_output.write('</svg>')


def draw_trna_domain_svg():
    """Visualize tRNA domains as a svg.

    """
    # using MT-TL1 as framework, so using genomic and pair coordinates here
    # first get information to draw lines between pairs
    matches_pairs = {}
    for row in csv.DictReader(open('required_files/svg_input/tRNA_model.tsv'), delimiter='\t'):
        if row["Pair_coordinate"]:
            matches_pairs[row["Pair_coordinate"]] = (row["x1"], row["y1"])
    
    # now build dictionary of domain per position
    matches_site = {}
    for row in csv.DictReader(
            open('output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt'), delimiter='\t'):
        if row["symbol"] == "MT-TL1":
            matches_site[row["tRNA_position"]] = row["tRNA_domain"]
    
    # use same colors as ggplot
    color_dict = {
        'Acceptorstem': '#F8766D', 'D-stem': '#CD9600', 'D-loop': '#7CAE00', 'Anticodonstem': '#00BE67',
        'Anticodonloop': '#00BFC4', 'Variableregion': '#00A9FF', 'T-stem': '#C77CFF', 'T-loop': '#FF61CC', '': 'white'}
    
    # now draw the svg
    with open('figure_scripts/extended_data_figures/FigureED6c_tRNA_domains.svg', 'w') as svg_output:
        svg_output.write('<svg height="350" width="350" xmlns="http://www.w3.org/2000/svg">\n\n')
        svg_output.write('\n<!-- MARKERS -->\n')
        
        for row in csv.DictReader(open('required_files/svg_input/tRNA_model.tsv'), delimiter='\t'):
            site = row["Site_number"]
            x1 = row["x1"]
            y1 = row["y1"]
            font = "arial"
            fontWeight = 'normal'
            fontSize = '7'
            lineColor = "#000000"
            strokeWidth = '1'
            radius = "5.5"  # circle radius
            textColor = "#000000"
            
            circleColor = color_dict[matches_site[site]]
            
            if ":a" in site:
                fontSize = '6'
                site = site.replace(':', "")
            
            circle = '<circle cx="' + str(x1) + '" cy="' + str(y1) + '" r="' + \
                     radius + '" style="fill: ' + circleColor + ';stroke:rgb(0,0,0);stroke-width:0.2" />'
            svg_output.write(circle)
            
            base_svg = '<text  text-anchor="middle" alignment-baseline="middle" x="' + x1 + '" y="' + y1 + \
                       '" style = "font-size: ' + fontSize + '; fill: ' + textColor + '; font-family: ' + font + \
                       '; font-weight: ' + fontWeight + ';" >' + site + '</text>'
            
            if row["Pair_coordinate"]:
                svg_output.write(base_svg)
                x1 = int(x1)
                x2 = int(matches_pairs[row["Genomic_coordinate"]][0])
                y1 = int(y1)
                y2 = int(matches_pairs[row["Genomic_coordinate"]][1])
                
                if x1 == x2:
                    x3 = x1
                    x4 = x2
                    middle_y = float((y1 + y2) / 2)
                    y3 = str((middle_y + 4))
                    y4 = str((middle_y - 4))
                elif y1 == y2:
                    middle_x = float((x1 + x2) / 2)
                    x3 = str((middle_x + 4))
                    x4 = str((middle_x - 4))
                    y3 = y1
                    y4 = y2
                else:
                    length_x = abs(float((x1 - x2)) / 3)
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
                
                line = '<line x1="' + str(x3) + '" y1="' + str(y3) + '" x2="' + str(x4) + '" y2="' + str(y4) + \
                       '" style="stroke: ' + lineColor + '; stroke-width:' + strokeWidth + '; stroke-linecap: round" />'
                svg_output.write(line)  # all as lines
            
            else:
                svg_output.write(base_svg)
        
        svg_output.write('</svg>')


if __name__ == "__main__":
    draw_trna_pos_svg()
    draw_trna_domain_svg()
