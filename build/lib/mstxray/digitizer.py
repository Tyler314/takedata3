try:
    import struck3302 as st
except:
    print('Error importing struck 3302 module.')

import time
import threading
import numpy as np


from remoteobject import serve_object, get_proxy

SRV_PORT = 18495
SRV_ADDR = 'xray.physics.wisc.edu'


class Server:
    def __init__(self):
        print('Initializing.')
        self.lock    = threading.RLock()
        self.handle  = st.init3150()
        self.cards   = None
        self.samples = 0
        print('Done.')


    def _get_avg(self, card, chan):
        with self.lock:
            self.arm(10000, False)
            self.trigger(False)
            time.sleep(0.1)
            data = self.get_data(card, chan, False)
            return data.mean()


    def _set_offset(self, card, chan, offset):
        with self.lock:
            st.set_dac_offset(self.handle, card, chan, int(offset))


    def set_cards(self, cardList):
        with self.lock:
            print('Setting the card list.')
            self.cards = cardList


    def zero(self, value=16384, delta=100):
        with self.lock:
            for card in self.cards:
                for chan in range(8):
                    print('Zeroing: card 0x{0:x} chan {1}.'.format(card, chan))
                    offset = value
                    self._set_offset(card, chan, offset)
                    avg = self._get_avg(card, chan)

                    while(np.abs(avg - value) > delta):
                        print('   Delta: {0}'.format(np.abs(avg - value)))
                        avg = self._get_avg(card, chan)
                        offset = offset + (value - avg)
                        self._set_offset(card, chan, offset)


    def arm(self, num_samples, verbose=True):
        with self.lock:

            # Make sure samples is a multiple of 4.
            num_samples -= num_samples % 4

            if verbose:
                print('Arming digitizer for {0} samples.'.format(num_samples))

            # The real number of samples digitized.
            self.samples = num_samples

            # Reset the cards.
            st.reset(self.handle, self.cards)

            # Set the acquisition settings.
            acq_list = [st.ACQ_DISABLE_AUTOSTART,
                        st.ACQ_DISABLE_MULTIEVENT,
                        st.ACQ_ENABLE_LEMO_START_STOP,
                        st.ACQ_DISABLE_INTERNAL_TRIGGER,
                        st.ACQ_CLOCK_LEMO]
            st.set_acquisition(self.handle, self.cards, acq_list)

            # Set the event configuration.
            st.set_event_config(self.handle, self.cards,
                                st.EVENT_ENABLE_SAMPLE_LENGTH_STOP)

            # Set the number of samples.
            st.set_samples(self.handle, self.cards, self.samples)

            # Set the trigger info.
            data = 0x5 + (0x10 << 8) + (0x10 << 16)
            st.set_trigger(self.handle, self.cards, data)

            # Threshold.
            data = 0x10000 - 0x100 + (0x1 << 24)
            st.set_trigger_thresh(self.handle, self.cards, data)

            # Set DAC values (ADC offsets) to mid range.
            # Clear first.
            st.set_dac_control_status(self.handle, self.cards, 3)

            # Arm the digitizer.
            st.arm(self.handle, self.cards)


    def trigger(self, verbose=True):
        with self.lock:
            if verbose:
                print('Triggered via software.')
            st.start_sampling(self.handle, self.cards)


    def triggered(self, verbose=True):
        with self.lock:
            if verbose:
                print('Checking if triggered.')
            return min(st.get_event_counter(self.handle, self.cards)) > 0


    def get_data(self, card, chan, verbose=True):
        with self.lock:
            if verbose:
                print('Reading: card 0x{0:x}, chan {1}.'.format(card, chan))
            arr = np.ndarray(self.samples, dtype=np.uint16)
            arr = np.ascontiguousarray(arr)
            st.read_channel_data(self.handle, card, chan, arr)
            return arr


def serve():
    serve_object(Server(), SRV_PORT)


def get_client():
    # Construct the Server proxy object.
    server = get_proxy(SRV_ADDR, SRV_PORT, Server)

    # Set the server card list. The list(set(...)) syntax removes
    # duplicates.
    import mstxray.physical as mp
    dets = mp.XRayDetector.get_all(in_use=True)
    cards = list(set([det.digiCard for det in dets]))
    server.set_cards(cards)

    return server
