import heapq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv('SiouxFalls_net.csv', delimiter="	")
nodes = [str(i) for i in range(1,25)]
links = {}
for i in range(len(df)):
    o = df.iloc[i,0]
    d = df.iloc[i,1]
    links[(str(o),str(d))] = {'capacity' : df.iloc[i,2], 'free_flow_time' : df.iloc[i,4]}

df = pd.read_csv('SiouxFalls_trips.csv', delimiter="	")


od_matrix = {}
for i in range(len(df)):
    o = str(df.iloc[i,0])
    d = str(df.iloc[i,1])
    if o == d:
        continue
    od_matrix[(o,d)] = df.iloc[i,2]
# Define OD matrix
init_od_matrix = {}
for node in nodes:
    for node1 in nodes:
        if node == node1:
            continue
        init_od_matrix[(node, node1)] = 0

for od in od_matrix:
    init_od_matrix[od] = od_matrix[od]
od_matrix = init_od_matrix

#link 초기화
for link in links:
    links[link]['travel_times'] = links[link]['free_flow_time']
    links[link]['marginal_travel_times'] = links[link]['free_flow_time']
    links[link]['flow'] = 0

# Link performance function (BPR)

def so_link_performance(flow, capacity, free_flow_time):
    alpha = 0.15
    beta = 4

    return free_flow_time * (1 + alpha * (flow / capacity) ** beta * (beta))


def link_performance(flow, capacity, free_flow_time):
    alpha = 0.15
    beta = 4
    return free_flow_time * (1 + alpha * (flow / capacity) ** beta )


# Initialize flow


def so_first_derivative(flow, capacity, free_flow_time):
    alpha = 0.15
    beta = 4
    return free_flow_time * (1 + alpha * beta * (flow / capacity) ** beta)

def so_second_derivative(flow, capacity, free_flow_time):
    alpha = 0.15
    beta = 4
    return free_flow_time * (1 + alpha * beta * (flow / capacity) ** beta)

def min_cap(links, path):
    m = float('inf')
    for od_pair in path:
        m = min(m, links[od_pair]['capacity'])
    return m
# Gradient projection method
def gradient_projection(links, od_matrix, max_iter=1000, tol=1e-4):
    Pi = {}
    for od_pair in od_matrix:
        Pi[od_pair] = dijkstra_shortest_path(nodes, links, od_pair[0], od_pair[1], so = True)
    #1.
    ttta=[]
    cola = 0
    while cola < 300:
        if cola%10 == 0:
            print(cola)
        for link in links:
            links[link]['flow'] = 0
        for od_pair in od_matrix:
            for pi in Pi[od_pair]:
                for link in pi[0]:
                    links[link]['flow'] += pi[1]
        for link in links:
            links[link]['travel_times'] = link_performance(links[link]['flow'], links[link]['capacity'], links[link]['free_flow_time'])
            links[link]['marginal_travel_times'] = so_link_performance(links[link]['flow'], links[link]['capacity'], links[link]['free_flow_time'])

        for od_pair in Pi:
            if od_matrix[od_pair] == 0:
                continue
            pistar = dijkstra_shortest_path(nodes, links, od_pair[0], od_pair[1], so=True)
            check = 1
            for pi in Pi[od_pair]:
                if pi[0] == pistar:
                    check = 0
                    break
            if check:
                Pi[od_pair].append([pistar,0])
            #3.
            if len(Pi[od_pair]) == 1:
                Pi[od_pair][0] = [Pi[od_pair][0][0], od_matrix[od_pair]]
            else:

                erase_list = []
                hsum = 0
                j = 0
                for i, pi in enumerate(Pi[od_pair]):
                    # print(pi[0], pistar)
                    if pi[0] == pistar:
                        j = i
                        continue
                    # c = path_cost(pi[0],links)
                    #
                    # cstar = path_cost(pistar, links)

                    xorpath = list(set(pi[0]) ^ set(pistar))
                    fd = 0

                    for od_pair1 in xorpath:
                        fd += so_first_derivative(links[od_pair1]['flow'], links[od_pair1]['capacity'], links[od_pair1]['free_flow_time'])
                    sd = 0
                    for od_pair1 in xorpath:
                        sd += so_second_derivative(links[od_pair1]['flow'], links[od_pair1]['capacity'], links[od_pair1]['free_flow_time'])
                    # print(Pi[od_pair][i][1])
                    Pi[od_pair][i][1] -= fd / sd * 0.95**cola
                    # print(c, cstar, sd)
                    Pi[od_pair][i][1] = max(0, Pi[od_pair][i][1])
                    Pi[od_pair][i][1] = min(min_cap(links,Pi[od_pair][i][0]), Pi[od_pair][i][1])

                    hsum += Pi[od_pair][i][1]
                    if Pi[od_pair][i][1] == 0:
                        erase_list.append(Pi[od_pair][i])
                Pi[od_pair][j][1] = od_matrix[od_pair]-hsum
                #4.
                # print(len(Pi[od_pair]))
                # print(Pi[od_pair])
                for er in erase_list:
                    Pi[od_pair].remove(er)
        ttt = 0
        for link in links:
            # print(link, links[link]['travel_times'], links[link]['flow'])
            ttt += links[link]['travel_times'] * links[link]['flow']
        ttta.append(ttt)



        # print(Pi[('1', '15')])
        cola += 1
    return ttta, links, Pi








    # for (origin, destination), demand in od_matrix.items():
    #     shortest_path = dijkstra_shortest_path(nodes, travel_times, origin, destination)
    #     for link in shortest_path:
    #         gradient[link] += travel_times[link] - travel_times[shortest_path[0]]

# Shortest path using Dijkstra's algorithm with a priority queue
import heapq


def dijkstra_shortest_path(nodes, links, origin, destination, so=False):
    pq = []
    heapq.heappush(pq, (0, origin))
    distances = {node: float('inf') for node in nodes}
    distances[origin] = 0
    previous_nodes = {node: None for node in nodes}
    path = []

    while pq:
        current_distance, current_node = heapq.heappop(pq)

        if current_node == destination:
            break

        for neighbor in nodes:
            if (current_node, neighbor) in links:
                if so:
                    distance = links[(current_node, neighbor)]['marginal_travel_times']
                else:
                    distance = links[(current_node, neighbor)]['travel_times']
                new_distance = current_distance + distance

                if new_distance < distances[neighbor]:
                    distances[neighbor] = new_distance
                    heapq.heappush(pq, (new_distance, neighbor))
                    previous_nodes[neighbor] = current_node

    current_node = destination
    while previous_nodes[current_node] is not None:
        path.insert(0, (previous_nodes[current_node], current_node))
        current_node = previous_nodes[current_node]

    return path, distances[destination]




# Run the gradient projection algorithm
def path_cost(path, links):
    cost = 0
    for link in path:
        cost = cost + links[link]['travel_times']
    return cost


ttta, links, Pi = gradient_projection(links,od_matrix)
for od_pair in Pi:
    if len(Pi[od_pair])>1:
        print(od_pair)
        print(Pi[od_pair])
        for p in Pi[od_pair]:
            ip = p[0]
            print(path_cost(ip,links))
plt.figure()
plt.plot(ttta)
plt.show()
print(ttta[len(ttta)-1])
# print(optimal_flow)
