# Script to run that allows us to run all of the major scripts of the program

def print_out():
    print "--MENU--"
    print "Enter 0 to generate models."
    print "Enter 1 to convert models to pdb files."
    print "Enter 2 to generate statistics file."
    print "Enter 3 to exit."      


while True:
    print_out()
    menu_input = raw_input()
    if int(menu_input.strip()) < 0 or int(menu_input.strip()) > 3:
        continue
    else:
        if int(menu_input.strip()) == 0:
            execfile("hamidMethod.py")      
        if int(menu_input.strip()) == 1:
            execfile("hamidGeneratePdb.py")
        if int(menu_input.strip()) == 2:
            execfile("hamidStat.py")
        if int(menu_input.strip()) == 3:
            exit()
