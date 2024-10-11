from utils import *
import json
import time
from datetime import datetime
from tc_commands import add_latency, add_latency_
import socket
from k8s import *
from kubernetes import client, config
import random
import threading
from krios_daemon import *
import numpy as np

satellites = read_tles("tles.txt")

config.load_kube_config()
api = client.CoreV1Api()
earth_radius = 6378.135
altitude = 550
theta = math.acos(earth_radius / (earth_radius + altitude))

ip_mapping = cloudlab_fetch_ip_mapping()

client_node = "jedi011"
ground_station = geodetic2cartesian(45.9174667, -119.2684488, 550000)
# leo_region_center = geodetic2cartesian(37.7749, -122.4194, 550000) #SF

# sat_mapping = {
#     # "jedi009" : 856,
#     # "jedi010" : 943,
#     "jedi053" : 1228,
#     "jedi058" : 1226,
#     # "sat1": 1226,
# }

def add_sat_labels(node):
    labels = get_node(api, node).metadata.labels
    val = str(sat_mapping[node])

    #need to do this to avoid the situation when rescheduler tries to access sat_id but it's deleted
    if labels.get("sat_id", None) != val:
        add_label_to_node(api, node, 'sat_id1', val)
        remove_label_from_node(api, node, 'sat_id')
        add_label_to_node(api, node, 'sat_id', val)
        remove_label_from_node(api, node, 'sat_id1')

    else:
        add_label_to_node(api, node, 'sat_id', val)

def update_sat_mapping(current_time, active_sat, node, passive_sat = None, passive_node = None):
    jd, fr = jday(current_time.year, current_time.month, current_time.day, current_time.hour, current_time.minute, current_time.second)

    sats = []
    maxval, mindist, optimal_sat = 0, 1000000000, -1
    for i, sat in enumerate(satellites):
        satellite = Satrec.twoline2rv(sat['line1'], sat['line2'])
        e, location, velocity = satellite.sgp4(jd, fr)
        val = (velocity[0] * (leo_region_center[0] - location[0])) + (velocity[1] * (leo_region_center[1] - location[1])) + (velocity[2] * (leo_region_center[2] - location[2]))
        dist = calculate_distance(location, leo_region_center)
        if dist < allowable_distance:
            if closest_node:
                if dist < mindist:
                    mindist = dist
                    optimal_sat = i
                    continue
            else:
                if maxval < val:
                    maxval = val
                    optimal_sat = i
                    continue
                
            if active_sat != i and passive_sat != i:
                sats.append((i, dist, val, location, velocity))

    visible_list = random.sample(range(0, len(sats)), min(3, len(sats)))
    i = 0
    for key in sat_mapping:
        if key == node or key == passive_node:
            continue

        if not random_node and i == 0 and optimal_sat != active_sat:
            sat_mapping[key] = optimal_sat
        else:
            sat_mapping[key] = sats[visible_list[i]][0]

        add_sat_labels(key)
        i = i + 1


pods = get_pods(api)
altitude = 550 # km
elevation_angle = np.radians(25) # radians
pod = pods[0]
radius = int(pod.metadata.labels.get('radius', 100))
app_location = parseLocation(pod.metadata.labels['leoregion'])
leo_region_center = geodetic2cartesian(app_location[0], app_location[1], 550000)
active_pod = pod.metadata.name
active_node = pod.spec.node_name
allowable_distance = get_allowable_distance(radius, altitude, elevation_angle)
worstcase = False
client_location = geodetic2cartesian(32.8801, -117.2340, 550000) # UCSD
# client_location = geodetic2cartesian(37.7749, -122.4194, 550000) # SF
orbit = -1
#if direction is true, satellites with ids lesser follow. If false, satellites with ids greater follow

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
host = ip_mapping[client_node]
port = 12444
print(host, port)
# s.connect((host,port))

expt_time = 6000
t1 = time.time()
random_node = False
closest_node = False
lookahead = True
msg = str(t1)
# s.send(msg.encode())
t1 = t1 + 4.5
tend = t1 + expt_time
start_time = datetime.fromtimestamp(t1)
jd, fr = jday(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)

sats = []
maxval, mindist, optimal_sat = 0, 1000000000, -1
for i, sat in enumerate(satellites):
    satellite = Satrec.twoline2rv(sat['line1'], sat['line2'])
    e, location, velocity = satellite.sgp4(jd, fr)
    val = (velocity[0] * (leo_region_center[0] - location[0])) + (velocity[1] * (leo_region_center[1] - location[1])) + (velocity[2] * (leo_region_center[2] - location[2]))
    dist = calculate_distance(location, leo_region_center)
    if dist < allowable_distance:
        if closest_node:
            if dist < mindist:
                mindist = dist
                active_sat = i
                continue
        else:
            if maxval < val:
                maxval = val
                active_sat = i
                if random_node:
                    sats.append((i, dist, val, location, velocity))
                continue

        sats.append((i, dist, val, location, velocity))

#print(maxval, active_sat)
print(sats)
visible_list = random.sample(range(0, len(sats)), min(4, len(sats)))
#print(len(sats), visible_list)
i = 0
if random_node:
    sat_mapping[active_node] = sats[visible_list[i]][0]
    i = i + 1
else:
    sat_mapping[active_node] = active_sat
add_sat_labels(active_node)

for key in sat_mapping:
    if key == active_node:
        continue

    sat_mapping[key] = sats[visible_list[i]][0]

    add_sat_labels(key)
    i = i + 1


print("all ready")

if time.time() < t1:
    time.sleep(t1 - time.time())

while True:
    if time.time() > tend:
        break

    while time.time() < t1:
        time.sleep(t1 - time.time())

    t1 = time.time()
    current_time = datetime.fromtimestamp(t1)
    #tc = threading.Thread(target = update_latency, args=(active_node,))
    #tc.start()

    print("main", t1, end = ' ')
    pods = get_pods(api)

    if len(pods) == 1:
        active_node = pods[0].spec.node_name
        update_sat_mapping(current_time, sat_mapping[active_node], active_node)
    elif len(pods) == 2:
        pod1, pod2 = pods[0], pods[1]
        if pods[0].metadata.creation_timestamp > pod2.metadata.creation_timestamp:
            pod1, pod2 = pods[1], pods[0]

        if is_pod_ready(api, pod2):
            active_node, active_pod, passive_node, passive_pod = pod2.spec.node_name, pod2.metadata.name, pod1.spec.node_name, pod1.metadata.name
            print("ready ", active_pod)
        else:
            active_node, active_pod, passive_node, passive_pod = pod1.spec.node_name, pod1.metadata.name, pod2.spec.node_name, pod2.metadata.name
            print("not ready ", passive_pod)

        update_sat_mapping(current_time, sat_mapping[active_node], active_node, sat_mapping[passive_node], passive_node)
    
    msg = active_node + "_" + json.dumps(sat_mapping)
    print(msg)
    s.send(msg.encode())
    
    t1 = t1 + 1

pods = get_pods(api)
for pod in pods:
    delete_pod(api, pod)
time.sleep(2)
msg = "exit"
print("exited")
s.send(msg.encode())
s.close()
re.join()
