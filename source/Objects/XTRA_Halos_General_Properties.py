class Halo_Status_Class():

   def __init__(self):
       # Halos information exist?
       self.HalosDataExist = False
       # AGNs information exist?
       self.AGNsDataExist = False
       # Solved for Lx, T, flux?
       self.LxTxSolved = False
       # Trasformed into XCat prefered coordinate?
       self.XCatPreferedCoordinate = False

   def update(self, Halo_data):
       self.HalosDataExist = True
       self.AGNsDataExist = False
       self.LxTxSolved = False
       self.XCatPreferedCoordinate = False
       
       


