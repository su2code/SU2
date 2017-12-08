import os, sys, shutil, copy, string

def convert_c2c(source_routine):

    source_routine += '_d'
    source_location = source_directory + '/' + source_routine + '.c'

    # check to see if converted file has already been created (if it exists, skip conversion)
    target_routine = source_routine

    target_location = target_directory + '/' + target_routine + '_debug.c'

    if os.path.isfile(target_location):
        print(target_routine + ' exists!')

    else:

        # open and read in source file
        source_file = open(source_location, 'r')

        source_file_data = source_file.readlines()

        # convert code and save to new_data
        new_data = []
        temp_data = []
        redefines = []

        # add new function definition
        for data_line in source_file_data:
            if data_line.strip().startswith('//SU2_C2CPP'):
                data_words = data_line.split()

                # redefine defines...
                if data_words[1] == 'REDEFINE':

                    redefines.append(data_words[2])
            #include code
            else:
                temp_data.append(data_line)

        # redefine defines...
        for data_line in temp_data:
            for redef in redefines:
                pos = redef.find('=')
                old = redef[:pos]
                new = redef[pos+1:]
                new = new[0] + new[2:]

                data_line = data_line.replace(old, new)

            new_data.append(data_line)

        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)


def convert_c2cpp(source_routine, overwrite=False):

    source_routine += '_d'
    source_directory = os.path.join(os.environ['SU2_HOME'],'SU2_CFD/src/c_routines_d')
    source_location = source_directory + '/' + source_routine + '.c'

    # check to see if converted file has already been created (if it exists, skip conversion)
    target_routine = source_routine

    target_directory = os.path.join(os.environ['SU2_HOME'],'SU2_CFD/src/cpp_routines_d')
    target_location = target_directory + '/' + target_routine + '.cpp'

    if os.path.isfile(target_location) and overwrite==False:
        print(target_routine + ' exists!')

    else:

        # open and read in source file
        source_file = open(source_location, 'r')

        source_file_data = source_file.readlines()

        # convert code and save to new_data
        new_data = []
        temp_data = []
        temp2_data = []
        redefines = []

        start = 0
        start_decl = 0
        counter = 0
        end_decl = 0

        for data_line in source_file_data:
            # check to open routine
            # check for open bracket
            if data_line.strip().endswith('{'):
                start = counter

            # check for cpp2c comment
            if data_line.strip().startswith('//SU2_C2CPP DIVIDER'):
                start_decl = start
                end_decl = counter

            counter = counter + 1


        # keep new vars from tapenade
        for data_line in source_file_data[start_decl+1:end_decl]:
            if data_line.strip().startswith('int') or data_line.strip().startswith('double'):
                data_words = data_line.split()
                if data_words[1].startswith('arg') or data_words[1].startswith('tmp') or data_words[1].startswith('ii'):
                    temp_data.append(data_line)

        # add new function definition
        for data_line in source_file_data[end_decl:]:
            if data_line.strip().startswith('//SU2_C2CPP'):
                data_words = data_line.split()

                if data_words[1] == 'INVARS': 
                    new_line = 'void ' + target_routine + '('
                    for word in data_words[2:]:
                        new_line += 'double ' + word + ', double ' + word + 'd, '

                elif data_words[1] == 'OUTVARS':
                    for word in data_words[2:]:
                        new_line += 'double ' + word + ', double ' + word + 'd, '
                    new_line = new_line[:-2]

                    new_line += ')\n{\n'
                    new_line += '//SU2_INSERT START\n'
                    temp_data.insert(0, new_line)

                # redefine defines...
                elif data_words[1] == 'REDEFINE':

                    redefines.append(data_words[2])
            #include rest of code
            else:
                temp_data.append(data_line)

        # redefine defines...
        for data_line in temp_data:
            for redef in redefines:
                pos = redef.find('=')
                old = redef[:pos]
                new = redef[pos+1:]
                new = new[0] + new[2:]

                data_line = data_line.replace(old, new)

            temp2_data.append(data_line)

        # add ending comment for insert

        end = 0
        counter = 0

        for data_line in temp2_data:
            if data_line.strip().endswith('}'):
                end = counter
            counter = counter + 1

        counter = 0
        for data_line in temp2_data:
            if counter == end:
                data_line = data_line.strip()
                data_line = data_line[:-1]
                data_line += '\n//SU2_INSERT END\n}'

            new_data.append(data_line)
            counter = counter + 1

        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)
        print('Wrote file: ' + target_location)

def convert_cpp2c(source_routine, source_location, overwrite=False):
    '''Convert a C++ code from within the SU2 source tree into a C code for automatic differentiation. The output will be located in $SU2_HOME/c_routines'''

    # Replace C++ :: notation with underscores
    target_routine = source_routine.replace('::', '__')

    # Location of output c file
    cfilename = target_routine + '.c'
    target_directory = os.path.join(os.environ['SU2_HOME'], 'SU2_CFD/src/c_routines')
    target_location  = os.path.join(target_directory,cfilename)

    # Do not process if file exists and overwrite is not set True
    if os.path.isfile(target_location) and overwrite==False:
        print(target_routine + ' exists!')

    else:

        # open and read in source file
        source_file = open(source_location, 'r')

        source_file_data = source_file.readlines()

        # convert code and save to new_data
        new_data = []
        routine_open = 0
        comment_open = 0
        inputs_open = 0
        declaration_open = 0
        subroutine_open = 0
        sub_open = 0
#        subroutine_list = []
        var_type = ''
        define_count = 0
        replace_line_open = 0

        for data_line in source_file_data:
            # check to open routine
            # check for cpp2c comment
            if data_line.strip().startswith('//SU2_CPP2C') and routine_open == 0 :
                # split into words
                data_words = data_line.split()

                # check for start of routine
                if data_words[1] == 'START':
                    if data_words[2] == source_routine:
                        # start declaration:
                        routine_open = 1
                        new_line = 'void ' + target_routine
                        #new_data.append(new_line)

            # if routine is open
            elif routine_open == 1:
                # check for cpp2c comment
                if data_line.strip().startswith('//SU2_CPP2C'):
                    # split into words
                    data_words = data_line.split()

                    # check to start/end commented section
                    if data_words[1] == 'COMMENT':
                        if data_words[2] == 'START':
                            comment_open = 1
                        elif data_words[2] == 'END':
                            comment_open = 0

                    # check for end of routine
                    elif data_words[1] == 'END':
                        if data_words[2] == source_routine:
                            # end routine
                            new_line = '}'
                            new_data.append(new_line)
                            routine_open = 0
                            break

                    # open/close input variable declaration
                    elif data_words[1] == 'CALL_LIST':
                        if (data_words[2] == 'START' and inputs_open == 0):
                            new_line += '('
                            extra_line = ''
                            inputs_open = 1
                        elif (data_words[2] == 'END' and inputs_open == 1):
                            new_line = new_line[:-2]
                            new_line += ')\r\n'
                            new_data.append(new_line)
                            inputs_open = 0

                            new_line = '{\r\n'
                            new_data.append(new_line)
                            divider_line = '\r\n//SU2_C2CPP DIVIDER\r\n' + extra_line
#                            new_data.append(new_line)
#                            new_data.append(extra_line)

#                            divider_line = extra_line

                    # get list of independent variables
                    elif (data_words[1] == 'INVARS' and inputs_open == 1):

                        extra_line += '//SU2_C2CPP INVARS'

                        var_type = 'double'

                        # declare inputs
                        for word in data_words[2:]:
                            new_line += var_type + ' ' + word + ', '
                            extra_line += ' ' + word

                        extra_line += '\r\n' 

                    # get list of dependent variables
                    elif (data_words[1] == 'OUTVARS' and inputs_open == 1):

                        extra_line += '//SU2_C2CPP OUTVARS'

                        var_type = 'double'

                        # declare outputs
                        for word in data_words[2:]:
                            new_line += var_type + ' ' + word + ', '
                            extra_line += ' ' + word

                        extra_line += '\r\n'

                    # get list of doubles/ints that need declaration (either input or internally declare)
                    elif data_words[1] == 'DECL_LIST':
                        if (data_words[2] == 'START' and declaration_open == 0):
                            declaration_open = 1
                        elif (data_words[2] == 'END' and declaration_open == 1):
                            declaration_open = 0
                            new_data.append(divider_line)

                    elif data_words[1] == 'VARS':
                        if data_words[2] == 'INT':
                            var_type = 'int'
                        elif data_words[2] == 'DOUBLE':
                            var_type = 'double'

                        if inputs_open == 1:
                            for word in data_words[3:]:
                                new_line += var_type + ' ' + word + ', '
                        elif declaration_open == 1:
                            if data_words[3] == 'SCALAR':
                                new_line = var_type + ' '
                                for word in data_words[4:]:
                                    new_line += word + ', '

                                new_line = new_line[:-2]
                                new_line += ';\r\n'
                                new_data.append(new_line)
                            elif data_words[3] == 'MATRIX':
                                dimensions = []
                                for word in data_words[4:]:
                                    if word.startswith('SIZE='):
                                        dimensions.append(word[5:])

                                new_line = ''

                                for word in data_words[4:]:
                                    if not word.startswith('SIZE='):
                                        new_line += var_type + ' '
                                        new_line += word

                                        for dim in dimensions:

                                            new_line += '[' + dim + ']'

                                        new_line +=  ';\r\n'

                                new_data.append(new_line)

                    # get #define constants
                    elif data_words[1] == 'DEFINE':
                        extra_line = ''
                        for word in data_words[2:]:
                            define_count = define_count + 1
                            extra_line += '#define ' + word + ' ' + str(1234000 + define_count) + '\r\n'
                            divider_line += '//SU2_C2CPP REDEFINE ' + str(1234000 + define_count) + '=' + word[0] + '_' + word[1:] + '\r\n' # inserts '_' in name to stop Tapenade replacing it with #define
                        new_data.insert(0, extra_line)

                    # get sub routine
                    elif data_words[1] == 'SUBROUTINE':
                        if data_words[2] == 'START' and subroutine_open == 0:
#                            new_line = data_words[3] + '('
                            subroutine_open = 1
                            sub_routine = data_words[3]
                        elif subroutine_open == 1:
                            if data_words[2] == 'END':
#                                new_line = new_line[:-2]
#                                new_line += ');\r\n'
#                                new_data.append(new_line)
                                new_sub_data = []
                                for sub_data_line in sub_data:
                                    var_count = 0
                                    for old_var in sub_vars_old:
                                        sub_data_line = sub_data_line.replace(old_var, sub_vars_new[var_count])

                                        var_count = var_count +1

                                    new_sub_data.append(sub_data_line)
                                for sub_data_line in new_sub_data:
                                    new_data.append(sub_data_line)
                                subroutine_open = 0
                            elif data_words[2] == 'LOCATION':
                                sub_location = source_directory + '/' + data_words[3]
                                sub_file = open(sub_location, 'r')
                                sub_file_data = sub_file.readlines()
                                sub_data = []
                                sub_vars_old = []

                                # now loop through to get sub routine code
                                for sub_data_line in sub_file_data:
                                    if sub_data_line.strip().startswith('//SU2_CPP2C') and sub_open == 0 :
                                        # split into words
                                        sub_data_words = sub_data_line.split()

                                        # check for start of subroutine
                                        if sub_data_words[1] == 'SUB':
                                            if sub_data_words[2] == 'START':
                                                if sub_data_words[3] == sub_routine:
                                                # start declaration:
                                                    sub_open = 1
                                    elif sub_open == 1:
                                        if sub_data_line.strip().startswith('//SU2_CPP2C'):
                                            # split into words
                                            sub_data_words = sub_data_line.split()

                                            if sub_data_words[1] == 'SUB':
                                                if sub_data_words[2] == 'END' and sub_data_words[3] == sub_routine:

                                                    sub_open = 0
                                                elif sub_data_words[2] == 'VARS':
                                                    for sub_var in sub_data_words[3:]:
                                                        sub_vars_old.append(sub_var)

                                        else:
                                            sub_data.append(sub_data_line)

                            elif data_words[2] == 'VARS':
                                sub_vars_new = []

                                for word in data_words[3:]:
                                    sub_vars_new.append(word)

                    # get new line
                    if data_words[1] == 'REPLACE_LINE':
                        if data_words[2] == 'START' and replace_line_open == 0:
                            replace_line_open = 1
                        elif data_words[2] == 'END' and replace_line_open == 1:
                            replace_line_open = 0
                        elif data_words[2] == 'NEW_LINE':
                            new_data_line = ''
                            for data_word in data_words[3:]:
                                new_data_line += data_word + ' '
                            new_data.append(new_data_line)

                    # get any required #includes (differentiable forms of non-differentiable standard functions)
#                    elif data_words[1] == 'INCLUDE':
#                        subroutine_known = 0
#                        for subroutine in subroutine_list:
#                            if (subroutine == data_words[2]):
#                                 subroutine_known = 1
#                        if (subroutine_known == 0):
#                            subroutine_list.append(data_words[2])

                # deal with non comment line
                elif comment_open == 0 and subroutine_open == 0 and replace_line_open == 0:
                    # cut and paste
                    new_data.append(data_line)

        # add in includes:
#        for subroutine in subroutine_list:
#            extra_line = '#include "' + subroutine + '.c"\r\n'
#            new_data.insert(0, extra_line)



        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)


def diff_routine(routine_name, invars, outvars, file_list):

    # Tapenade command
    tapenade = 'tapenade '

    # Tapenade mode
    mode = '-tangent '

    # Routine
    root = '-root ' + routine_name + ' '

    # Variables
    independents = '-vars "' + ' '.join(invars) + '"' + ' '
    dependents = '-outvars "' + ' '.join(outvars) + '"' + ' '

    # Output directory
    output_directory = '-outputdirectory ' + os.path.join(os.environ['SU2_HOME'],'SU2_CFD/src/c_routines_d') + ' '
    # Files
    files = ' '.join(file_list)

    # Send command
    tapenade_call = tapenade + mode + root + independents + dependents + \
        output_directory + files
    print('Differentiating ' + routine_name)
    print('Command: $ ' + tapenade_call)
    os.system(tapenade_call)


def insert_math(source_routine):

    source_routine += '_d'
    source_location = source_directory + '/' + source_routine + '.c'

    # open and read in source and target files
    source_file = open(source_location, 'r')

    source_file_data = source_file.readlines()

    header_file = open(header_location, 'r')

    header_file_data = header_file.readlines()

    target_file = open(target_location, 'r')

    target_file_data = target_file.readlines()

    new_data = []
    temp_data = []
    temp_data2 = []
    routine_open = 0
    comment_open = 0


    # get function to add in
    for data_line in source_file_data:
        if data_line.strip().startswith('/*'):
            comment_open = 1
        elif data_line.strip().endswith('*/'):
            comment_open = 0
        elif comment_open == 0:
            #include code
            temp_data.append(data_line)

    for data_line in temp_data:
        if not data_line.strip().startswith('//'):
            temp_data2 = data_line
            break

    # add in header
    for data_line in header_file_data:
        if data_line.strip().startswith('//SU2_DIFF'):

            new_data.append(data_line)

            data_words = data_line.split()

            if data_words[1] == 'START' and routine_open == 0:
                if data_words[2] == source_routine:
                    routine_open = 1

                    new_data.append('\n')

                    new_data_line = temp_data2.replace('{',';') 

                    new_data.append(new_data_line)

                    new_data.append('\n')

            elif data_words[1] == 'END' and routine_open == 1:
                if data_words[2] == source_routine:
                    routine_open = 0

        elif routine_open == 0:
            new_data.append(data_line)

        # open and write header file
        header_file = open(header_location, 'w')

        header_file.writelines(new_data)

    new_data = []

    # add in function
    for data_line in target_file_data:
        if data_line.strip().startswith('//SU2_DIFF'):

            new_data.append(data_line)

            data_words = data_line.split()

            if data_words[1] == 'START' and routine_open == 0:
                if data_words[2] == source_routine:
                    routine_open = 1

                    new_data.append('\n')

                    for new_data_line in temp_data:
                        new_data.append(new_data_line)

                    new_data.append('\n')

            elif data_words[1] == 'END' and routine_open == 1:
                if data_words[2] == source_routine:
                    routine_open = 0

        elif routine_open == 0:
            new_data.append(data_line)

        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)

def insert_routine(source_routine, target_routine, target_file):

    source_routine += '_d'
    source_location = source_directory + '/' + source_routine + '.cpp'

    target_location = target_directory + '/' + target_file


    # open and read in source and target files
    source_file = open(source_location, 'r')

    source_file_data = source_file.readlines()

    target_file = open(target_location, 'r')

    target_file_data = target_file.readlines()

    new_data = []
    temp_data = []
    routine_open = 0

    # get function to add in
    for data_line in source_file_data:
        if data_line.strip().startswith('//SU2_INSERT'):
            data_words = data_line.split()

            if data_words[1] == 'START' and routine_open == 0:
                routine_open = 1
            elif data_words[1] == 'END' and routine_open == 1:
                routine_open = 0
        elif routine_open == 1:      
            #include code
            temp_data.append(data_line)

    routine_open = 0

    # add in function
    for data_line in target_file_data:
        if data_line.startswith('//SU2_DIFF'):

            new_data.append(data_line)

            data_words = data_line.split()

            if data_words[1] == 'START' and routine_open == 0:
                if data_words[2] == target_routine:
                    routine_open = 1

                    new_data.append('\n')

                    for new_data_line in temp_data:
                        new_data.append(new_data_line)

                    new_data.append('\n')

            elif data_words[1] == 'END' and routine_open == 1:
                if data_words[2] == target_routine:
                    routine_open = 0

        elif routine_open == 0:
            new_data.append(data_line)

        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)
