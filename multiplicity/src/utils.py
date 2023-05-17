def replace_elements(input_string, pars):
        # Replace elements in the input string with brackets
        for element in pars:
            input_string = input_string.replace(element, f'[{element}]')
        
        if "[" not in input_string:
            raise ValueError("ERROR: Unable to find variables in custom bin string")
            # If no brackets are found in the modified string, raise an error
        
        return input_string