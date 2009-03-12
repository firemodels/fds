"""
Here is a silly example of its usage:

import progress
import time
import random

total = 1000
p = progress.ProgressMeter(total=total)

while total > 0:
    cnt = random.randint(1, 25)
    p.update(cnt)
    total -= cnt
    time.sleep(random.random())


Here is an example of its output:

[------------------------->                                   ] 41%  821.2/sec


2006-02-20 Denis Barmenkov: ANSI codes replaced by Backspace (0x08) characters
"""
import time, sys, math

class ProgressMeter(object):
    #ESC = chr(27)
    def __init__(self, **kw):
        # What time do we start tracking our progress from?
        self.timestamp = kw.get('timestamp', time.time())
        # What kind of unit are we tracking?
        self.unit = str(kw.get('unit', ''))
        # Number of units to process
        self.total = int(kw.get('total', 100))
        # Number of units already processed
        self.count = int(kw.get('count', 0))
        # Refresh rate in seconds
        self.rate_refresh = float(kw.get('rate_refresh', .5))
        # Number of ticks in meter
        self.meter_ticks = int(kw.get('ticks', 60))
        self.meter_division = float(self.total) / self.meter_ticks
        self.meter_value = int(self.count / self.meter_division)
        self.last_update = None
        self.rate_history_idx = 0
        self.rate_history_len = 10
        self.rate_history = [None] * self.rate_history_len
        self.rate_current = 0.0
        self.last_refresh = 0
        self.prev_meter_len = 0

    def update(self, count, **kw):
        now = time.time()
        # Caclulate rate of progress
        rate = 0.0
        # Add count to Total
        self.count += count
        self.count = min(self.count, self.total)
        if self.last_update:
            delta = now - float(self.last_update)
            if delta:
                rate = count / delta
            else:
                rate = count
            self.rate_history[self.rate_history_idx] = rate
            self.rate_history_idx += 1
            self.rate_history_idx %= self.rate_history_len
            cnt = 0
            total = 0.0
            # Average rate history
            for rate in self.rate_history:
                if rate == None:
                    continue
                cnt += 1
                total += rate
            rate = total / cnt
        self.rate_current = rate
        self.last_update = now
        # Device Total by meter division
        value = int(self.count / self.meter_division)
        if value > self.meter_value:
            self.meter_value = value
        if self.last_refresh:
            if (now - self.last_refresh) > self.rate_refresh or \
                (self.count >= self.total):
                    self.refresh()
        else:
            self.refresh()

    def get_meter(self, **kw):
        bar = '-' * self.meter_value
        pad = ' ' * (self.meter_ticks - self.meter_value)
        perc = (float(self.count) / self.total) * 100
        return '[%s>%s] %d%%  %.1f/sec' % (bar, pad, perc, self.rate_current)

    def refresh(self, **kw):
        # Clear line and return cursor to start-of-line
        sys.stdout.write(' ' * self.prev_meter_len + '\x08' * self.prev_meter_len)
        # Get meter text
        meter_text = self.get_meter(**kw)
        # Write meter and return cursor to start-of-line
        sys.stdout.write(meter_text + '\x08'*len(meter_text))
        self.prev_meter_len = len(meter_text)

        # Are we finished?
        if self.count >= self.total:
            sys.stdout.write('\n')
        sys.stdout.flush()
        # Timestamp
        self.last_refresh = time.time()
