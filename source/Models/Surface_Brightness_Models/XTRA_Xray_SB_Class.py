class Surface_Brightness_Class():

   def __init__(self, option=1):

      import sys 
      if (option==1):
         from simpleBetaProfile import Surface_Brightness_Model
      elif (option==2): 
         from simpleDoubleBetaProfile import Surface_Brightness_Model
      elif (option==3): 
         from simpleCCandNCCProfile import Surface_Brightness_Model
      else:
         print "ERROR: Input option does not exist!"
         print "Please change XML file and try again."
         raw_input("Press enter to exit! ")
         sys.exit(2)
      
      self.F2SB_Class = Surface_Brightness_Model()

