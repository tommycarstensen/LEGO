import xml.etree.ElementTree as ET
import random

def get_bounds(bricks_section):
    x_min, x_max = float('inf'), float('-inf')
    z_min, z_max = float('inf'), float('-inf')

    for Brick in bricks_section.findall('Brick'):
        for Part in Brick.findall('Part'):
            Bone = Part.find('Bone')
            if Bone is not None:
                transformation = list(map(float, Bone.get('transformation').split(',')))
                x, z = transformation[9], transformation[11]
                x_min = min(x_min, x)
                x_max = max(x_max, x)
                z_min = min(z_min, z)
                z_max = max(z_max, z)

    return x_min, x_max, z_min, z_max

def calculate_percentage(start_value, end_value, position, start_position, end_position):
    if position <= start_position:
        return start_value
    elif position >= end_position:
        return end_value
    else:
        return start_value + (end_value - start_value) * (position - start_position) / (end_position - start_position)

def apply_color_transitions(input_lxfml, output_lxfml, bottom_id, middle_id, top_id):
    tree = ET.parse(input_lxfml)
    root = tree.getroot()

    bricks_section = root.find('Bricks')
    if bricks_section is None:
        print("No Bricks section found in the LXFML file.")
        return

    # Get dynamic boundaries for horizontal gradients
    x_min, x_max, z_min, z_max = get_bounds(bricks_section)

    for Brick in list(bricks_section.findall('Brick')):
        for Part in list(Brick.findall('Part')):
            Bone = Part.find('Bone')
            if Bone is None:
                continue

            transformation = list(map(float, Bone.get('transformation').split(',')))
            x_coordinate = transformation[9]
            y_coordinate = transformation[10]
            z_coordinate = transformation[11]
            layer = int(round(y_coordinate / 0.96, 0))

            if layer < 15:
                bricks_section.remove(Brick)
                break

            if 15 <= layer <= 59:
                # Vertical gradient (based on layer)
                bottom_percentage = calculate_percentage(1.0, 0.0, layer, 15, 44)
                middle_percentage = 0.0
                if 15 <= layer <= 29:
                    middle_percentage = calculate_percentage(0.0, 0.5, layer, 15, 29)
                elif 30 <= layer <= 44:
                    middle_percentage = 0.5
                elif 45 <= layer <= 59:
                    middle_percentage = calculate_percentage(0.5, 0.0, layer, 45, 59)
                top_percentage = calculate_percentage(0.0, 1.0, layer, 30, 44)

                # Horizontal gradient (x and z influence)
                x_gradient = calculate_percentage(1.0, 0.0, x_coordinate, x_min, x_max)
                z_gradient = calculate_percentage(1.0, 0.0, z_coordinate, z_min, z_max)

                # Combine gradients (e.g., average or multiply, depending on desired effect)
                combined_bottom_percentage = bottom_percentage * x_gradient * z_gradient
                combined_middle_percentage = middle_percentage * x_gradient * z_gradient
                combined_top_percentage = top_percentage * x_gradient * z_gradient

                # Normalize percentages to avoid too much influence from one component
                total_percentage = (combined_bottom_percentage + combined_middle_percentage + combined_top_percentage)
                if total_percentage == 0:
                    total_percentage = 1  # Prevent division by zero

                combined_bottom_percentage /= total_percentage
                combined_middle_percentage /= total_percentage
                combined_top_percentage /= total_percentage

                rand_value = random.random()
                if rand_value <= combined_bottom_percentage:
                    materialID = bottom_id
                elif rand_value <= combined_bottom_percentage + combined_middle_percentage:
                    materialID = middle_id
                else:
                    materialID = top_id

                Part.set('materials', f"{materialID},0")

    tree.write(output_lxfml, encoding='UTF-8', xml_declaration=True)
    print(f"New LXFML file created: {output_lxfml}")

if __name__ == "__main__":
    input_lxfml = 'input.lxfml'

    # For Ice Hands
    apply_color_transitions(
        input_lxfml,
        'output_ice.lxfml',
        top_id="42",  # ice: trans-light blue
        middle_id="40",  # ice: trans-clear
        bottom_id="43"  # ice: trans-dark blue
    )

    # For Fire Hands
    apply_color_transitions(
        input_lxfml,
        'output_fire.lxfml',
        bottom_id="44",  # fire: trans-yellow
        middle_id="182",  # fire: trans-orange
        top_id="41"  # fire: trans-red
    )
