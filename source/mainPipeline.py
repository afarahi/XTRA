import sys
import time

def makingXRayCatalogs():

    import random
    from Objects import Halos_Class

    random.seed(int(sys.argv[1]))

    ##################################################################
    #@@@@ STEP 1 [Initializing]
    #@@@@
    sTime = time.time()
    Halos = Halos_Class()
    eTime = time.time()-sTime
    print "Eplased time to initialize the code : ", eTime, " s"

    ##################################################################
    #@@@@ STEP 2 [Reading halos]
    #@@@@
    sTime = time.time()
    Halos.addHalosCatalog(fname=sys.argv[2])
    eTime = time.time()-sTime
    print "Eplased time to read halo catalogs : ", eTime, " s"

    ##################################################################
    #@@@@ STEP 3 [Solving for LX and Tx]
    #@@@@
    sTime = time.time()
    Halos.SolveLxTxFlux()
    eTime = time.time()-sTime
    print "Eplased time to assign X-ray observables : ", eTime, " s"

    ##################################################################
    #@@@@ STEP 4 [Saving clusters]
    #@@@@
    sTime = time.time()
    Halos.SaveHalosData()
    eTime = time.time()-sTime
    print "Eplased time to save outputs : ", eTime, " s"


def makingXRayRealization():

    from Objects import Map_Halos_Class
    from Objects import Map_Class

    ##################################################################
    #@@@@ STEP 1 [Initializing and reading clusters catalogs]
    #@@@@
    sTime = time.time()
    clusters = Map_Halos_Class()
    clusters.update(fname=sys.argv[2])
    eTime = time.time()-sTime
    print "Eplased time to initialize and read halo catalog: ",eTime, " s"

    ##################################################################
    #@@@@ STEP 2 [Reading halos]
    #@@@@
    sTime = time.time()
    Map = Map_Class()
    Map.update()
    eTime = time.time()-sTime
    print "Eplased time to initialize the map object: ", eTime, " s"

    ##################################################################
    #@@@@ STEP 3 [Adding clusters to the surface brightness map]
    #@@@@
    sTime = time.time()
    Map.addClusters2Map(clusters)
    eTime = time.time()-sTime
    print "Eplased time to add clusters to the surface brightness map: ", eTime, " s"
    
    if 0: 
      # WKB method: 
      sTime = time.time()
      Map.drawHaloBounds(clusters) 
      eTime = time.time()-sTime
      print "Elapsed time for drawing:",eTime," s"

    ##################################################################
    #@@@@ STEP 4 [Saving maps]
    #@@@@
    sTime = time.time()
    if 0: 
      Map.saveMapPic()
    else: 
      Map.saveMapPic(draw_halos=True,Halos=clusters) # WKB
    Map.saveMapFits()
    eTime = time.time()-sTime
    print "Eplased time to save map: ", eTime, " s"


def makingEventMap():
    
    from Objects import Event_Map_Class

    ##################################################################
    #@@@@ STEP 1 [Initializing the event map]
    #@@@@ 
    sTime = time.time()
    EMap = Event_Map_Class()
    eTime = time.time()-sTime
    print "Eplased time to initialize event map: ", eTime, " s"

    ##################################################################
    #@@@@ STEP 2 [Generate one tiled realization]
    #@@@@
    sTime = time.time()
    EMap.run_tiled_map(fname=sys.argv[2])
    eTime = time.time()-sTime
    print "Eplased time to initialize event map: ", eTime, " s"





