import os,sys

def get_channel_data( chanfile, data_dir="." ):
    """ collect channel information from files channels folder """
    
    channeldata = { "channame":[],
                    "channum":[],
                    "cpstate":[],
                    "flavor":[],
                    "num_target_factor":[],
                    "numchans":0}

    chanfilename = data_dir+"/channels/channels_"+chanfile+".dat";
    if not os.path.exists(chanfilename):
        raise ValueError( "could not find channel file: ",chanfilename)
    
    chans = open(chanfilename,'r')
    for l in chans.readlines():
        l = l.strip()
        info = l.split()
        channeldata["channame"].append( info[0] )
        channeldata["channum"].append( int(info[1]) )
        channeldata["cpstate"].append( info[2] )
        channeldata["flavor"].append(  info[3] )
        channeldata["num_target_factor"].append( float(info[4]) )
        channeldata["numchans"] += 1
    chans.close();

    print("Number of channels: ",channeldata["numchans"])
    return channeldata
