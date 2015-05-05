#from genie import *
#from genie_epics_api import *
import genie_python.genie as gen
import genie_python.genie_epics_api as geapi

# genie_epics_api is needed for set_sample_par amongst other things

def dosans(size):
    print "dosans"
    setuplarmor_normal()

def dotrans():
    print "dotrans"

def move():
    print "moving sample"

    
def do_hist_sans():
    print "do_hist_sans"

def setuplarmor_scanning():
    gen.change(nperiods=1)
    print "setup larmor scanning"
    gen.change_start()
    gen.change_tables(spectra="C:\Instrument\Settings\Tables\spectra_scanning_80.dat")
    gen.change_tables(wiring="C:\Instrument\Settings\Tables\wiring.dat")
    gen.change_tcb(low=5.0,high=100000.0,step=100.0,trange=1)
    gen.change_finish()

def setuplarmor_nr():
    gen.change(nperiods=1)
    print "setup larmor nr"
    gen.change_start()
    gen.change_tables(spectra="C:\Instrument\Settings\Tables\spectra_nrscanning.dat")
    gen.change_tables(wiring="C:\Instrument\Settings\Tables\wiring.dat")
    gen.change_tcb(low=5.0,high=100000.0,step=100.0,trange=1)
    gen.change_finish()

def setuplarmor_nrscanning():
    gen.change(nperiods=1)
    print "setup larmor nrscanning"
    gen.change_start()
    #change_tables(spectra="C:\Instrument\Settings\Tables\spectra_nrscanning.dat")
    gen.change_tables(spectra="U:\Users\Masks\spectra_scanning_auto.dat")
    gen.change_tables(wiring="C:\Instrument\Settings\Tables\wiring.dat")
    gen.change_tcb(low=5.0,high=100000.0,step=100.0,trange=1)
    gen.change_finish()

def setuplarmor_normal():
    print "setup larmor normal"
    gen.change(nperiods=1)
    gen.change_start()
    gen.change_tables(spectra="C:\Instrument\Settings\Tables\spectra_1To1.dat")
    gen.change_tables(wiring="C:\Instrument\Settings\Tables\wiring.dat")
    gen.change_tcb(low=5.0,high=100000.0,step=100.0,trange=1)
    gen.change_finish()

def setuplarmor_polarised():
    print "setup larmor polarised"
    gen.change(nperiods=1)
    gen.change_start()
    gen.change_tables(spectra="C:\Instrument\Settings\Tables\spectra_1To1.dat")
    gen.change_tables(wiring="C:\Instrument\Settings\Tables\wiring.dat")
    gen.change_tcb(low=5.0,high=100000.0,step=100.0,trange=1)
    gen.change_finish()

def setuplarmor_bsalignment():
    gen.change(nperiods=1)
    print "setup larmor beamstop alignment"
    gen.change_start()
    gen.change_tables(spectra="C:\Instrument\Settings\Tables\spectra_1To1.dat")
    gen.change_tables(wiring="C:\Instrument\Settings\Tables\wiring.dat")
    gen.change_tcb(low=1000.0,high=100000.0,step=99000.0,trange=1)
    gen.change_finish()

def setuplarmor_monitorsonly():
    gen.change(nperiods=1)
    print "setup larmor monitors only"
    gen.change_start()
    gen.change_tables(spectra="C:\Instrument\Settings\Tables\spectra_phase1.dat")
    gen.change_tables(wiring="C:\Instrument\Settings\Tables\wiring.dat")
    gen.change_tcb(low=1000.0,high=100000.0,step=100.0,trange=1)
    gen.change_finish()

def setuplarmor_event():
    gen.change(nperiods=1)
    print "setup larmor for event mode"
    gen.change_start()
    gen.change_tables(spectra="C:\Instrument\Settings\Tables\spectra_1To1.dat")
    gen.change_tables(wiring="C:\Instrument\Settings\Tables\wiring_event.dat")
    gen.change_tcb(low=5.0,high=100000.0,step=100.0,trange=1)
    gen.change_finish()

def setuplarmor_4periods():
    gen.change(nperiods=1)
    print "setup larmor for 4 Period mode"
    gen.change_start()
    gen.change_tables(spectra="C:\Instrument\Settings\Tables\spectra_4To1.dat")
    gen.change_tables(wiring="C:\Instrument\Settings\Tables\wiring.dat")
    gen.change_tcb(low=5.0,high=100000.0,step=100.0,trange=1)
    gen.change_finish()

def setuplarmor_quiet():
    print "setup larmor quiet"

def FOMin():
    cset(pol_trans=110,pol_arc=-1.6)

def ShortPolariserin():
    cset(pol_trans=9,pol_arc=-1.6)

def LongPolariserin():
    cset(pol_trans=209,pol_arc=-1.6)
    
def homecoarsejaws():
    print "Homing Coarse Jaws"
    gen.cset(cjhgap=40,cjvgap=40)
    gen.waitfor_move()
    # home north and west
    gen.set_pv("IN:LARMOR:MOT:JAWS1:JN:MTR.HOMR","1")
    gen.set_pv("IN:LARMOR:MOT:JAWS1:JW:MTR.HOMR","1")
    gen.waitfor_move()
    gen.set_pv("IN:LARMOR:MOT:JAWS1:JN:MTR.VAL","20")
    gen.set_pv("IN:LARMOR:MOT:JAWS1:JW:MTR.VAL","20")
    # home south and east
    gen.set_pv("IN:LARMOR:MOT:JAWS1:JS:MTR.HOMR","1")
    gen.set_pv("IN:LARMOR:MOT:JAWS1:JE:MTR.HOMR","1")
    gen.waitfor_move()
    gen.set_pv("IN:LARMOR:MOT:JAWS1:JS:MTR.VAL","20")
    gen.set_pv("IN:LARMOR:MOT:JAWS1:JE:MTR.VAL","20")
    

def homea1():
    print "Homing a1"

def homes1():
    print "Homing s1"

def homes2():
    print "Homing s2"

def detectoronoff(onoff=0,delay=1):
    if(onoff==1):
        gen.set_pv("IN:LARMOR:CAEN:hv0:0:8:pwonoff","On")
        gen.set_pv("IN:LARMOR:CAEN:hv0:0:9:pwonoff","On")
        gen.set_pv("IN:LARMOR:CAEN:hv0:0:10:pwonoff","On")
        gen.set_pv("IN:LARMOR:CAEN:hv0:0:11:pwonoff","On")
    else:
        gen.set_pv("IN:LARMOR:CAEN:hv0:0:8:pwonoff","Off")
        gen.set_pv("IN:LARMOR:CAEN:hv0:0:9:pwonoff","Off")
        gen.set_pv("IN:LARMOR:CAEN:hv0:0:10:pwonoff","Off")
        gen.set_pv("IN:LARMOR:CAEN:hv0:0:11:pwonoff","Off")
    # wait for 3 minutes for ramp up 60s to ranmp down
    if(delay==1):
        if(onoff==1):
            print "Waiting For Detector To Power Up (180s)"
            sleep(180)
        else:
            print "Waiting For Detector To Power Down (60s)"
            sleep(60)
            

def movebench(angle=0.0,delaydet=1):
    print "Turning Detector Off"
    detectoronoff(onoff=0,delay=delaydet)
    a1=0
    a1+=gen.get_pv("IN:LARMOR:CAEN:hv0:0:8:status")
    a1+=gen.get_pv("IN:LARMOR:CAEN:hv0:0:9:status")
    a1+=gen.get_pv("IN:LARMOR:CAEN:hv0:0:10:status")
    a1+=gen.get_pv("IN:LARMOR:CAEN:hv0:0:11:status")
    if(a1>0):
        print "The detector is not turned off"
        print "Not attempting Move"
        return
    else:
        print "The detector is off"
       
    if(angle >= 0.0):
        gen.cset(benchlift=1)
        print "Lifting Bench (20s)"
        sleep(20)
        a1=gen.get_pv("IN:LARMOR:BENCH:STATUS")

        if(a1==1):
            print "Rotating Bench"
            gen.cset(bench_rot=angle)
            waitfor_move()
            print "Lowering Bench (20s)"
            gen.cset(benchlift=0)
            sleep(20)
        else:
            print "Bench failed to lift"
            print "Move not attempted"
    #turn the detector back on
    print "Turning Detector Back on"
    detectoronoff(onoff=1,delay=delaydet)

def rotatebench(angle=0.0,delaydet=1):
    a1=0
    a1+=gen.get_pv("IN:LARMOR:CAEN:hv0:0:8:status")
    a1+=gen.get_pv("IN:LARMOR:CAEN:hv0:0:9:status")
    a1+=gen.get_pv("IN:LARMOR:CAEN:hv0:0:10:status")
    a1+=gen.get_pv("IN:LARMOR:CAEN:hv0:0:11:status")
    if(a1>0):
        print "The detector is not turned off"
        print "Not attempting Move"
        return
    else:
        print "The detector is off"
       
    if(angle >= -0.5):
        gen.cset(benchlift=1)
        print "Lifting Bench (20s)"
        sleep(20)
        a1=gen.get_pv("IN:LARMOR:BENCH:STATUS")

        if(a1==1):
            print "Rotating Bench"
            gen.cset(bench_rot=angle)
            waitfor_move()
            print "Lowering Bench (20s)"
            gen.cset(benchlift=0)
            sleep(20)
        else:
            print "Bench failed to lift"
            print "Move not attempted"

