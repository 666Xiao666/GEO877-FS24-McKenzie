import heapq
import itertools
import matplotlib.pyplot as plt
import math

# Define the Point class to represent a point in 2D space
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

# Define the Event class to represent events in the Fortune's algorithm
class Event:
    def __init__(self, x, p, a):
        self.x = x
        self.p = p
        self.a = a
        self.valid = True

# Define the Arc class to represent an arc in the beach line
class Arc:
    def __init__(self, p, a=None, b=None):
        self.p = p
        self.pprev = a
        self.pnext = b
        self.e = None
        self.s0 = None
        self.s1 = None

# Define the Segment class to represent the edges of the Voronoi diagram
class Segment:
    def __init__(self, p):
        self.start = p
        self.end = None
        self.done = False

    # Method to mark the segment as finished by setting its end point
    def finish(self, p):
        if self.done: return
        self.end = p
        self.done = True        

# Define the PriorityQueue class for managing the priority queue of events and points
class PriorityQueue:
    def __init__(self):
        self.pq = []
        self.entry_finder = {}
        self.counter = itertools.count()

    # Method to add an item to the priority queue
    def push(self, item):
        if item in self.entry_finder: return
        count = next(self.counter)
        entry = [item.x, count, item]
        self.entry_finder[item] = entry
        heapq.heappush(self.pq, entry)

    # Method to remove an item from the priority queue
    def remove_entry(self, item):
        entry = self.entry_finder.pop(item)
        entry[-1] = 'Removed'

    # Method to pop the item with the highest priority (lowest value) from the queue
    def pop(self):
        while self.pq:
            priority, count, item = heapq.heappop(self.pq)
            if item != 'Removed':
                del self.entry_finder[item]
                return item
        raise KeyError('pop from an empty priority queue')

    # Method to get the item with the highest priority without removing it from the queue
    def top(self):
        while self.pq:
            priority, count, item = heapq.heappop(self.pq)
            if item != 'Removed':
                del self.entry_finder[item]
                self.push(item)
                return item
        raise KeyError('top from an empty priority queue')

    # Method to check if the priority queue is empty
    def empty(self):
        return not self.pq

# Define the Voronoi class to implement Fortune's algorithm for generating Voronoi diagrams
class Voronoi:
    def __init__(self, points, bbox):
        self.output = []    # List to store the output segments
        self.arc = None     # Pointer to the current arc in the beach line

        self.points = PriorityQueue()    # Priority queue for points
        self.event = PriorityQueue()     # Priority queue for circle events

        self.x0, self.x1, self.y0, self.y1 = bbox  # Bounding box for the diagram

        # Initialize the points in the priority queue
        for pts in points:
            point = Point(pts[0], pts[1])
            self.points.push(point)

    # Method to process the events and points to generate the Voronoi diagram
    def process(self):
        # Process all points and events
        while not self.points.empty():
            if not self.event.empty() and (self.event.top().x <= self.points.top().x):
                self.process_event()
            else:
                self.process_point()

        # Process remaining events
        while not self.event.empty():
            self.process_event()

        # Finish the edges of the Voronoi diagram
        self.finish_edges()

    # Method to process a site event (point)
    def process_point(self):
        p = self.points.pop()
        self.arc_insert(p)

    # Method to process a circle event
    def process_event(self):
        e = self.event.pop()
        if e.valid:
            s = Segment(e.p)
            self.output.append(s)
            a = e.a
            if a.pprev is not None:
                a.pprev.pnext = a.pnext
                a.pprev.s1 = s
            if a.pnext is not None:
                a.pnext.pprev = a.pprev
                a.pnext.s0 = s

            if a.s0 is not None: a.s0.finish(e.p)
            if a.s1 is not None: a.s1.finish(e.p)

            if a.pprev is not None: self.check_circle_event(a.pprev, e.x)
            if a.pnext is not None: self.check_circle_event(a.pnext, e.x)

    # Method to insert an arc into the beach line
    def arc_insert(self, p):
        if self.arc is None:
            self.arc = Arc(p)
        else:
            i = self.arc
            while i is not None:
                flag, z = self.intersect(p, i)
                if flag:
                    flag, zz = self.intersect(p, i.pnext)
                    if (i.pnext is not None) and (not flag):
                        i.pnext.pprev = Arc(i.p, i, i.pnext)
                        i.pnext = i.pnext.pprev
                    else:
                        i.pnext = Arc(i.p, i)
                    i.pnext.s1 = i.s1

                    i.pnext.pprev = Arc(p, i, i.pnext)
                    i.pnext = i.pnext.pprev

                    i = i.pnext

                    seg = Segment(z)
                    self.output.append(seg)
                    i.pprev.s1 = i.s0 = seg

                    seg = Segment(z)
                    self.output.append(seg)
                    i.pnext.s0 = i.s1 = seg

                    self.check_circle_event(i, p.x)
                    self.check_circle_event(i.pprev, p.x)
                    self.check_circle_event(i.pnext, p.x)

                    return
                i = i.pnext

            i = self.arc
            while i.pnext is not None:
                i = i.pnext
            i.pnext = Arc(p, i)

            x = self.x0
            y = (i.pnext.p.y + i.p.y) / 2.0
            start = Point(x, y)

            seg = Segment(start)
            i.s1 = i.pnext.s0 = seg
            self.output.append(seg)

    # Method to check and add a circle event for an arc
    def check_circle_event(self, i, x0):
        if (i.e is not None) and (i.e.x != self.x0):
            i.e.valid = False
        i.e = None

        if (i.pprev is None) or (i.pnext is None): return

        flag, x, o = self.circle(i.pprev.p, i.p, i.pnext.p)
        if flag and (x > self.x0):
            i.e = Event(x, o, i)
            self.event.push(i.e)

    # Method to compute the circumcircle of three points
    def circle(self, a, b, c):
        if ((b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y)) > 0: return False, None, None

        A = b.x - a.x
        B = b.y - a.y
        C = c.x - a.x
        D = c.y - a.y
        E = A*(a.x + b.x) + B*(a.y + b.y)
        F = C*(a.x + c.x) + D*(a.y + c.y)
        G = 2*(A*(c.y - b.y) - B*(c.x - b.x))

        if (G == 0): return False, None, None

        ox = 1.0 * (D*E - B*F) / G
        oy = 1.0 * (A*F - C*E) / G

        x = ox + math.sqrt((a.x-ox)**2 + (a.y-oy)**2)
        o = Point(ox, oy)
           
        return True, x, o
        
    # Method to check if a new point intersects an arc
    def intersect(self, p, i):
        if (i is None): return False, None
        if (i.p.x == p.x): return False, None

        a = 0.0
        b = 0.0

        if i.pprev is not None:
            a = (self.intersection(i.pprev.p, i.p, 1.0*p.x)).y
        if i.pnext is not None:
            b = (self.intersection(i.p, i.pnext.p, 1.0*p.x)).y

        if (((i.pprev is None) or (a <= p.y)) and ((i.pnext is None) or (p.y <= b))):
            py = p.y
            px = 1.0 * ((i.p.x)**2 + (i.p.y-py)**2 - p.x**2) / (2*i.p.x - 2*p.x)
            res = Point(px, py)
            return True, res
        return False, None

    # Method to find the intersection point of two parabolas
    def intersection(self, p0, p1, l):
        p = p0
        if (p0.x == p1.x):
            py = (p0.y + p1.y) / 2.0
        elif (p1.x == l):
            py = p1.y
        elif (p0.x == l):
            py = p0.y
            p = p1
        else:
            z0 = 2.0 * (p0.x - l)
            z1 = 2.0 * (p1.x - l)

            a = 1.0/z0 - 1.0/z1
            b = -2.0 * (p0.y/z0 - p1.y/z1)
            c = 1.0 * (p0.y**2 + p0.x**2 - l**2) / z0 - 1.0 * (p1.y**2 + p1.x**2 - l**2) / z1

            py = 1.0 * (-b-math.sqrt(b*b - 4*a*c)) / (2*a)
            
        px = 1.0 * (p.x**2 + (p.y-py)**2 - l**2) / (2*p.x-2*l)
        res = Point(px, py)
        return res

    # Method to finish the edges of the Voronoi diagram
    def finish_edges(self):
        l = self.x1 + (self.x1 - self.x0) + (self.y1 - self.y0)
        i = self.arc
        while i.pnext is not None:
            if i.s1 is not None:
                p = self.intersection(i.p, i.pnext.p, l*2.0)
                i.s1.finish(p)
            i = i.pnext

    # Method to get the output segments of the Voronoi diagram
    def get_output(self):
        res = []
        for o in self.output:
            p0 = o.start
            p1 = o.end
            if p1 is not None:
                res.append((p0.x, p0.y, p1.x, p1.y))
        return res

# Function to plot the Voronoi diagram using matplotlib
def plot_voronoi(points, segments):
    fig, ax = plt.subplots()
    for (x1, y1, x2, y2) in segments:
        ax.plot([x1, x2], [y1, y2], 'r')
    px, py = zip(*points)
    ax.plot(px, py, 'bo')
    ax.set_xlim(548000, 575000)
    ax.set_ylim(250000, 270000)
    plt.show()

# Define the input points and bounding box
points = [(556854.0927575239, 258756.71976257954), (552961.2590527163, 260673.33194784168), 
          (557751.7340146792, 253393.05717845634), (558969.2728328079, 259127.1139168013), 
          (557067.1893649272, 254941.66919528414), (557860.1872208333, 257161.78357722983), 
          (557393.4972634824, 258584.8512827782), (556977.183444837, 255725.65688890964), 
          (557041.5507356384, 256757.29464616627), (556195.7373974773, 258775.9757858105), 
          (556555.9063358027, 253831.93935211096), (558932.496566676, 255747.44426909648), 
          (562783.6789938582, 254844.56230746955), (557728.3067704869, 254221.14371476788), 
          (556782.881363944, 255534.73743720353), (553758.6597355691, 259155.8399370797), 
          (557455.8177479993, 255837.85328022018), (557351.2761499344, 261484.2865470387), 
          (558129.9393348523, 256418.46992108598), (558261.2771765835, 258282.52602208406), 
          (556946.733090108, 257939.94360529166), (557728.8561302526, 255140.799285179), 
          (557513.5748439579, 257686.31171812583), (556596.5562943101, 255959.66234214418), 
          (557175.7621433879, 255284.6547875153), (557044.9879401936, 257340.36001075804), 
          (556362.4076066122, 257340.81948037166), (556719.7898231349, 256721.58446211927), 
          (557245.270139444, 255631.3134805858), (561328.7589288311, 255740.70473248605), 
          (557503.2180973361, 255132.7755585173), (557003.6959896815, 256413.03932245448), 
          (560067.1815628994, 259921.51731126942), (565085.4064372301, 266583.1113133952), 
          (551593.5896089211, 261398.56842957158), (556366.1894671313, 255116.41707014944), 
          (555118.0147069613, 257298.6738503268), (556549.2518765404, 258407.98592053074), 
          (553554.5629041132, 258220.67992737982), (565048.3493186166, 266125.14927044325), 
          (552473.4208245791, 261990.3818048099), (557311.2275252102, 254300.5387907531), 
          (555909.9493602022, 258084.2024959987), (556215.8818559591, 262362.16619821545), 
          (557269.5120540475, 260071.84092636872), (554550.5224943471, 257636.4710407704), 
          (556425.9321898739, 261477.17045246437), (555026.1043057269, 258284.97925114818), 
          (557204.260873609, 256353.88684566412), (556579.6773157837, 260562.78456048295), 
          (560599.8367114541, 256325.63292417116), (557426.3939656046, 254806.199898554), 
          (564611.9772507473, 265243.84609774407), (556596.8324934612, 256303.90641559474), 
          (558453.1247341277, 255420.88668588363), (564685.3029210464, 265573.7993946895), 
          (552096.5062934769, 259787.72563414834), (558155.2673504568, 255966.34614766855), 
          (558380.0523823582, 253124.95972373057), (553943.006695772, 258166.27314490918), 
          (556756.4159311918, 259272.6731852591), (553169.7328056679, 257864.91276576743), 
          (559112.2793790943, 254433.58318478148), (557213.4304828062, 256052.38145246077), 
          (560085.5090686884, 261267.6279687956), (556888.5897497272, 258207.08978696167), 
          (557303.2578624435, 253904.44197674748), (556643.8764576514, 259996.83898436185), 
          (557143.3901693798, 256572.70199671015), (562435.778566557, 255079.7014344195), 
          (556139.2999334992, 259194.89988034312), (563636.9314088896, 255028.04045851808), 
          (555622.5809521555, 260587.54299806524), (558218.7652923979, 253678.3866825141), 
          (559964.993043698, 255133.63231500145), (558159.8285596772, 255562.90595101006), 
          (558184.5140792694, 255114.55333461147), (558448.0981708481, 255928.68718734756), 
          (553196.5018601543, 259505.92447243258), (562069.2490681551, 256635.0584571343), 
          (565222.2673860164, 265940.5054450361), (565022.3763443012, 265098.4946719734), 
          (565372.2671699915, 264773.86068887543), (564665.2355105373, 263967.4405574212)]
bbox = (548000, 250000, 575000, 270000)

# Create an instance of the Voronoi class and process the points
fortune = Voronoi(points, bbox)
fortune.process()
output_segments = fortune.get_output()

print(output_segments)
# Plot the resulting Voronoi diagram
plot_voronoi(points, output_segments)
