import json

with open('XTRA_info.json') as fp:
    _info = json.load(fp)

__author__ = _info['author']
__version__ = _info['version']
__revision__ = _info['revision']
__revnumber__ = int(__revision__.strip('$')[9:])
__date__ = _info['date']

del _info

__logo__ = " ##############################################################\n"
__logo__ += " ##                                                 V.%s  ##\n"%__version__
__logo__ += " ##  XX      XX     XXXXXXXXXXXX   XXXXXXX         xXXx      ##\n"
__logo__ += " ##   XX    XX           XX        XX    XX       XXxxXX     ##\n"
__logo__ += " ##    XX  XX            XX        XX     XX     XX    XX    ##\n"
__logo__ += " ##     XxxX             XX        XXXXXXXX     XX      XX   ##\n"
__logo__ += " ##     xXXx   XXXX      XX        XXXxx        XX      XX   ##\n"
__logo__ += " ##     XxxX             XX        XX  XX       XX XXXX XX   ##\n"
__logo__ += " ##    XX  XX            XX        XX   XX      XX      XX   ##\n"
__logo__ += " ##   XX    XX           XX        XX    XX     XX      XX   ##\n"
__logo__ += " ##  XX      XX          XX        XX     XX    XX      XX   ##\n"
__logo__ += " ##                                            C ARYA FARAHI ##\n"
__logo__ += " ##############################################################\n"

