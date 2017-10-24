import numpy as np

def getData4(derp):
        channelData = []
        HEADER_LENGTH = 3
        FILE_HEADER = 16

        # Read in the .dat file
        try:
            data = open(str(derp) + '.dat', 'r');
        except:
            print str(derp) + ".dat not found"
            sys.exit(1)

        # Store every 2 bytes (16 bits) into an element of the array 'f'
        f = np.fromfile(data, dtype=np.uint16)
        # "start" is the starting point of the current set of data points being evaluated
        start = 10
        # "offset" is the starting point of the header associated with the current set of
        # data points being evaluated.
        offset = 0
	# Set a counter for the channels
	current_channel = 0
	print "len(f) = " + str(len(f))
        # Iterate through the entire set of channels and their corresponding data points.
        while len(f) >= offset + 11:
            print 'CHANNEL {}: '.format(current_channel),
            print hex(f[offset] & 0xffff)[0:-1] + ', ' + hex(f[offset + 1] & 0xffff)[0:-1]
            
            # Initialize a list to hold the data points of the desired channel
            data_points = []
            # Grab the lower 2 bytes of the number indicating how many samples there are for the channel
            lowerBits = f[offset + 8]
            # Grab the upper 2 bytes of this number
            upperBits = f[offset + 9]
            # Create the number
            sample_length = 2 * (((upperBits & 0xfff) << 16) + lowerBits)
            print "Event Length = {}".format(sample_length)
	    # Check if valid channel, if not then exit
            if len(f) < (sample_length + start):
                sys.exit()
            # store the data points of this channel
            for i in range(start, sample_length + start):
                data_points.append(f[i])
            # Store the data points as a numpy array, contiguous in memory
            channelData.append(np.ascontiguousarray(np.asarray(data_points)))
 
            # Update the starting point of the data, the header, and the channel number
            start += (sample_length + 10)
            offset += (start - 10)
            current_channel += 1

        return channelData
