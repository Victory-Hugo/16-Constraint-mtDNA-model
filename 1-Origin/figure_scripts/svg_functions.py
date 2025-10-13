from typing import IO, List

# building upon code found online to create manual gradient
# see https://bsouthga.dev/posts/color-gradients-with-python


def hex_to_RGB(hex: str):
	"""Turns hex color into rgb format, ie #FFFFFF -> [255,255,255]
	
	:param hex: color in hex format
	:return: color in rgb format, tuple of integers
	"""
	# Pass 16 to the integer function for change of base
	return [int(hex[i:i + 2], 16) for i in range(1, 6, 2)]


def RGB_to_hex(rgb: List[int]):
	"""Turns rgb color into hex format, ie [255,255,255] -> #FFFFFF
	
	:param rgb: color in rgb format, tuple of integers
	:return: color in hex format
	"""
	# Components need to be integers for hex to make sense
	rgb = [int(x) for x in rgb]
	return "#" + "".join(["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in rgb])


def color_dict(gradient: List[List[int]]):
	"""Takes in a list of rgb sublists and returns dictionary of colors in rgb and hex form for use in plotting
	
	:param gradient: list of rgb colors
	:return: dictionary with hex, r, g, b as keys, and ordered list of color values for each key
	"""
	return {
		"hex": [RGB_to_hex(rgb) for rgb in gradient], "r": [rgb[0] for rgb in gradient],
		"g": [rgb[1] for rgb in gradient], "b": [rgb[2] for rgb in gradient]}


def linear_gradient(start_hex: str, finish_hex: str, n: int):
	"""Returns a gradient list of (n) colors between two hex colors
	
	:param start_hex: start of gradient color in format (ie "#FFFFFF")
	:param finish_hex: same format as start_hex, for end of gradient
	:param n: number of colors to include in gradient list
	:return: dictionary with hex, r, g, b as keys, and ordered list of n color values for each key
	"""
	# Starting and ending colors in RGB form
	s = hex_to_RGB(start_hex)
	f = hex_to_RGB(finish_hex)
	# Initialize a list of the output colors with the starting color
	RGB_list = [s]
	# Calculate a color at each evenly spaced value of t from 1 to n
	for t in range(1, n):
		# Interpolate RGB vector for color at the current value of t
		curr_vector = [
			int(s[j] + (float(t) / (n - 1)) * (f[j] - s[j]))
			for j in range(3)
		]
		# Add it to our list of output colors
		RGB_list.append(curr_vector)
	
	return color_dict(RGB_list)


def return_color(value: float):
	"""Determine the fill color to represent constraint in svg plots, white to red gradient on log scale

	:param value: the obs:exp ratio upper CI (OEUF)
	:return: the fill color to use, in hex format
	"""
	
	# blended gradient
	def merge(dict1, dict2):
		for k, v in dict2.items():
			if k in dict1:
				dict1[k].extend(v)
			else:
				dict1[k] = v
		return dict1
	
	dict = merge(
		linear_gradient(start_hex="#e41a1c", finish_hex="#fbe7e7", n=75),
		linear_gradient(start_hex="#fbe7e7", finish_hex="#fdfafa", n=31))
	
	if value >= 1.05:
		fill_color = dict['hex'][-1]
	else:
		fill_color = dict['hex'][round(value * 100)]
		
	return fill_color


def color_legend(svg_output: IO, x: int, y1: int, width: int, height: int, font_size: int):
	"""Print a legend of the color gradient used for visualizing constraint

	:param svg_output: svg to append the legend to
	:param x: x coordinate of the legend
	:param y1: y coordinate for the top of the legend
	:param width: width of the rectangle color legend
	:param height: height of the rectangle color legend
	:param font_size: font size in legend
	"""
	# apply + 5% offset here to place 1.0 correctly
	svg_output.write(
		'<defs>\n \
			<linearGradient id="LegendGradient" gradientTransform="rotate(90)">\n \
				<stop offset="0%" stop-color="' + "#fdfafa" + '" />\n \
				<stop offset="25%" stop-color="' + "#fbe7e7" + '" />\n \
				<stop offset="105%" stop-color="' + "#e41a1c" + '" />\n \
			</linearGradient>\n \
		</defs>\n')
	
	svg_output.write(
		'<rect x="' + str(x) + '" y="' + str(y1) + '" width= "' + str(width) + '" height= "' + str(height) +
		'" fill="' + 'url(#LegendGradient)' + '"/>')
	svg_output.write(
		'<text text-anchor="start" x="' + str(x) + '" y="' + str(y1 - 20) + '" style = "font-size: ' +
		str(font_size + 3) + '; fill: #000000; font-family: arial; font-weight: normal;" > '
		'OEUF legend </text>')
	
	# add legend ticks and labels
	# using y1 + y, where y = height - (value * height)
	# y = height - ((value / max_value) * height) if use gradient beyond 1
	svg_output.write(
		'<text text-anchor="start" x="' + str(x + width) + '" y="' + str(y1 + height + 2) + '" style = "font-size: ' +
		str(font_size + 2) + '; fill: #000000; font-family: arial; font-weight: normal;" > '
		'- 0.0 </text>')
	svg_output.write(
		'<text text-anchor="start" x="' + str(x + width) + '" y="' + str(y1 + (height - ((1/1.05) * height))) + '" style = "font-size: ' +
		str(font_size + 2) + '; fill: #000000; font-family: arial; font-weight: normal;" > '
		'- 1.0 </text>')
	svg_output.write(
		'<text text-anchor="start" x="' + str(x + width) + '" y="' + str(y1 + height - ((0.5 / 1.05) * height)) + '" style = "font-size: ' +
		str(font_size + 2) + '; fill: #000000; font-family: arial; font-weight: normal;" > - ' + str(0.5) + '</text>')
