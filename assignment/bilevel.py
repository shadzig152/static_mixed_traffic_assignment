import heapq
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import copy


df = pd.read_csv('SiouxFalls_net.csv', delimiter="	")
nodes = []
links = {}
for i in range(len(df)):
    nodes.append(str(df.iloc[i,0]))
    o = df.iloc[i,0]
    d = df.iloc[i,1]
    links[(str(o),str(d))] = {'capacity' : df.iloc[i,2], 'free_flow_time' : df.iloc[i,4]}
nodes = list(set(nodes))
df = pd.read_csv('SiouxFalls_trips.csv', delimiter="	")

od_matrix = {}
for i in range(len(df)):
    o = str(df.iloc[i,0])
    d = str(df.iloc[i,1])
    if o == d:
        continue
    if df.iloc[i,2] == 0:
        continue
    od_matrix[(o,d)] = df.iloc[i,2]
# Define OD matrix


# od_matrix[('1', '20')] *= 10




# Link performance function (BPR)
def link_performance(flow, capacity, free_flow_time):
    alpha = 0.15
    beta = 4
    return free_flow_time * (1 + alpha * (flow / capacity) ** beta)


def so_link_performance(flow, capacity, free_flow_time):
    alpha = 0.15
    beta = 4

    return free_flow_time * (1 + alpha * (flow / capacity) ** beta * (beta))


# Initialize flow

def path_cost(path, links):
    cost = 0
    for link in path:
        cost = cost + links[link]['travel_times']
    return cost
def second_derivative(link, links):
    alpha = 0.15
    beta = 4
    return links[link]['free_flow_time']*alpha*beta*(links[link]['flow']*1.0/links[link]['capacity'])**(beta-1)/links[link]['capacity']

def so_first_derivative(flow, capacity, free_flow_time):
    alpha = 0.15
    beta = 4
    return free_flow_time * (1 + alpha * beta * (flow / capacity) ** beta)

def so_second_derivative(flow, capacity, free_flow_time):
    alpha = 0.15
    beta = 4
    return free_flow_time * alpha * beta * beta / capacity * (flow/capacity)**(beta-1)


# Gradient projection method
def gradient_projection(links, od_matrix, e):
    Pi_ue = {}
    Pi_so = {}
    so_od_matrix = {}
    ue_od_matrix = {}
    # link 초기화
    for link in links:
        links[link]['travel_times'] = links[link]['free_flow_time']
        links[link]['marginal_travel_times'] = links[link]['free_flow_time']
        links[link]['flow'] = 0


    #1.
    for od in od_matrix:
        Pi_ue[od] = []
        Pi_so[od] = []


    ################# ue 100% #######################


    ue_od_matrix = copy.deepcopy(od_matrix)

    cola = 0
    while cola < 300:

        ###################################
        for od_pair in Pi_ue:
            ###################ue flow 분배########################
            if ue_od_matrix[od_pair] != 0:

                pistar, cost = dijkstra_shortest_path(nodes, links, od_pair[0], od_pair[1], so=False)
                check = 1
                for pi in Pi_ue[od_pair]:
                    if pi[0] == pistar:
                        check = 0
                        break
                if check:
                    Pi_ue[od_pair].append([pistar, 0])
                # 3.
                if len(Pi_ue[od_pair]) == 1:
                    Pi_ue[od_pair][0] = [Pi_ue[od_pair][0][0], ue_od_matrix[od_pair]]
                else:
                    erase_list = []
                    hsum = 0
                    j = 0
                    for i, pi in enumerate(Pi_ue[od_pair]):
                        if pi[0] == pistar:
                            j = i
                            break
                    for i, pi in enumerate(Pi_ue[od_pair]):
                        # print(pi[0], pistar)
                        if pi[0] == pistar:
                            continue
                        c = path_cost(pi[0], links)

                        cstar = path_cost(pistar, links)

                        xorpath = list(set(pi[0]) ^ set(pistar))
                        sd = 0
                        for od_pair1 in xorpath:
                            sd += second_derivative(od_pair1, links)
                        # print(Pi[od_pair][i][1])
                        Pi_ue[od_pair][i][1] -= (c - cstar) / sd * 0.13 * 0.98 ** cola
                        # print(c, cstar, sd)
                        Pi_ue[od_pair][i][1] = max(0, Pi_ue[od_pair][i][1])
                        hsum += Pi_ue[od_pair][i][1]
                        if Pi_ue[od_pair][i][1] == 0:
                            erase_list.append(Pi_ue[od_pair][i])
                    Pi_ue[od_pair][j][1] = ue_od_matrix[od_pair] - hsum
                    # 4.
                    # print(len(Pi[od_pair]))
                    # print(Pi[od_pair])
                    for er in erase_list:
                        Pi_ue[od_pair].remove(er)
            ################################################################

        ######### link 업데이트#############
        for link in links:
            links[link]['flow'] = 0
        for od_pair in od_matrix:
            for pi in Pi_ue[od_pair]:
                for link in pi[0]:
                    links[link]['flow'] += pi[1]
        for link in links:
            links[link]['travel_times'] = link_performance(links[link]['flow'], links[link]['capacity'],
                                                           links[link]['free_flow_time'])
            links[link]['marginal_travel_times'] = so_link_performance(links[link]['flow'],
                                                                       links[link]['capacity'],
                                                                               links[link]['free_flow_time'])



        cola += 1


    ####### od 나누기 #######
    for od in od_matrix:

        so_od_matrix[od] = od_matrix[od]*e
        ue_od_matrix[od] = od_matrix[od]*(1-e)

    ########################

    Pi_so = copy.deepcopy(Pi_ue)

    for od in Pi_ue:
        for j, pi in enumerate(Pi_ue[od]):
            Pi_so[od][j][1] = Pi_ue[od][j][1]*e
            Pi_ue[od][j][1] *= (1-e)




    step = 0.002
    # step = 0.002
    cola = 0
    while cola < 300:

        ###################################
        for od_pair in Pi_ue:
            ###################so flow 분배########################
            if so_od_matrix[od_pair] != 0:
                pistar, cost = dijkstra_shortest_path(nodes, links, od_pair[0], od_pair[1], so=True)
                check = 1
                for pi in Pi_so[od_pair]:
                    if pi[0] == pistar:
                        check = 0
                        break
                if check:
                    Pi_so[od_pair].append([pistar, 0])
                # 3.
                if len(Pi_so[od_pair]) == 1:
                    Pi_so[od_pair][0] = [Pi_so[od_pair][0][0], so_od_matrix[od_pair]]
                else:
                    erase_list = []
                    hsum = 0
                    j = 0
                    for i, pi in enumerate(Pi_so[od_pair]):
                        if pi[0] == pistar:
                            j = i
                            break
                    for i, pi in enumerate(Pi_so[od_pair]):
                        # print(pi[0], pistar)
                        if pi[0] == pistar:
                            continue

                        # for od in pi[0]:
                        #     fd1 += so_first_derivative(links[od]['flow'], links[od]['capacity'], links[od]['free_flow_time'])
                        # for od in pistar:
                        #     fd2 += so_first_derivative(links[od]['flow'], links[od]['capacity'], links[od]['free_flow_time'])
                        fd = 0
                        xorpath = list(set(pi[0]) ^ set(pistar))

                        for od_pair1 in xorpath:
                            fd += so_first_derivative(links[od_pair1]['flow'], links[od_pair1]['capacity'],
                                                      links[od_pair1]['free_flow_time'])


                        sd = 0
                        for od_pair1 in xorpath:
                            sd += so_second_derivative(links[od_pair1]['flow'], links[od_pair1]['capacity'], links[od_pair1]['free_flow_time'])
                        # print(Pi[od_pair][i][1])
                        Pi_so[od_pair][i][1] -= fd / sd * step * 0.99 ** cola
                        # print(c, cstar, sd)
                        Pi_so[od_pair][i][1] = max(0, Pi_so[od_pair][i][1])
                        hsum += Pi_so[od_pair][i][1]
                        if Pi_so[od_pair][i][1] == 0:
                            erase_list.append(Pi_so[od_pair][i])
                    Pi_so[od_pair][j][1] = so_od_matrix[od_pair] - hsum
                    # 4.
                    # print(len(Pi[od_pair]))
                    # print(Pi[od_pair])
                    for er in erase_list:
                        Pi_so[od_pair].remove(er)
            ###############################################################

        ######### link 업데이트#############

        for link in links:
            links[link]['flow'] = 0
        for od_pair in od_matrix:
            for pi in Pi_ue[od_pair]:
                for link in pi[0]:
                    links[link]['flow'] += pi[1]
            for pi in Pi_so[od_pair]:
                for link in pi[0]:
                    links[link]['flow'] += pi[1]

        for link in links:
            links[link]['travel_times'] = link_performance(links[link]['flow'], links[link]['capacity'],
                                                           links[link]['free_flow_time'])
            links[link]['marginal_travel_times'] = so_link_performance(links[link]['flow'],
                                                                       links[link]['capacity'],
                                                                       links[link]['free_flow_time'])


        cola +=1

    cola = 0
    while cola < 300:

        for od_pair in Pi_ue:
            ###################ue flow 분배########################
            if ue_od_matrix[od_pair] != 0:

                pistar, cost = dijkstra_shortest_path(nodes, links, od_pair[0], od_pair[1], so=False)
                check = 1
                for pi in Pi_ue[od_pair]:
                    if pi[0] == pistar:
                        check = 0
                        break
                if check:
                    Pi_ue[od_pair].append([pistar, 0])
                # 3.
                if len(Pi_ue[od_pair]) == 1:
                    Pi_ue[od_pair][0] = [Pi_ue[od_pair][0][0], ue_od_matrix[od_pair]]
                else:
                    erase_list = []
                    hsum = 0
                    j = 0
                    for i, pi in enumerate(Pi_ue[od_pair]):
                        if pi[0] == pistar:
                            j = i
                            break
                    for i, pi in enumerate(Pi_ue[od_pair]):
                        # print(pi[0], pistar)
                        if pi[0] == pistar:
                            continue
                        c = path_cost(pi[0], links)

                        cstar = path_cost(pistar, links)

                        xorpath = list(set(pi[0]) ^ set(pistar))
                        sd = 0
                        for od_pair1 in xorpath:
                            sd += second_derivative(od_pair1, links)
                        # print(Pi[od_pair][i][1])
                        Pi_ue[od_pair][i][1] -= (c - cstar) / sd * step * 0.99 ** cola
                        # print(c, cstar, sd)
                        Pi_ue[od_pair][i][1] = max(0, Pi_ue[od_pair][i][1])
                        hsum += Pi_ue[od_pair][i][1]
                        if Pi_ue[od_pair][i][1] == 0:
                            erase_list.append(Pi_ue[od_pair][i])
                    Pi_ue[od_pair][j][1] = ue_od_matrix[od_pair] - hsum
                    # 4.
                    # print(len(Pi[od_pair]))
                    # print(Pi[od_pair])
                    for er in erase_list:
                        Pi_ue[od_pair].remove(er)
            ################################################################

        ######### link 업데이트#############
        for link in links:
            links[link]['flow'] = 0
        for od_pair in od_matrix:
            for pi in Pi_ue[od_pair]:
                for link in pi[0]:
                    links[link]['flow'] += pi[1]
            for pi in Pi_so[od_pair]:
                for link in pi[0]:
                    links[link]['flow'] += pi[1]
        for link in links:
            links[link]['travel_times'] = link_performance(links[link]['flow'], links[link]['capacity'],
                                                           links[link]['free_flow_time'])
            links[link]['marginal_travel_times'] = so_link_performance(links[link]['flow'],
                                                                       links[link]['capacity'],
                                                                       links[link]['free_flow_time'])
        ###################################
        cola += 1





    # print(Pi[('1', '15')])
    ttt=0
    for link in links:
        # print(link, links[link]['travel_times'], links[link]['flow'])
        ttt += links[link]['travel_times'] * links[link]['flow']

    print(str(e) + ':' + str(ttt))
    # ttt = 0
    # for link in links:
    #     # print(link, links[link]['travel_times'], links[link]['flow'])
    #     ttt += links[link]['travel_times'] * links[link]['flow']
    return ttt, Pi_ue, Pi_so, copy.deepcopy(links)








    # for (origin, destination), demand in od_matrix.items():
    #     shortest_path = dijkstra_shortest_path(nodes, travel_times, origin, destination)
    #     for link in shortest_path:
    #         gradient[link] += travel_times[link] - travel_times[shortest_path[0]]

# Shortest path using Dijkstra's algorithm with a priority queue
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


elog = []

plt.figure()
data = {}
for i in range(11):
    e = i/10
    ttt, Pi_ue, Pi_so, output_links = gradient_projection(copy.deepcopy(links), od_matrix, e)

    for link in links:
        links[link]['flow'] = 0
    for od_pair in od_matrix:
        for pi in Pi_ue[od_pair]:
            for link in pi[0]:
                links[link]['flow'] += pi[1]
        for pi in Pi_so[od_pair]:
            for link in pi[0]:
                links[link]['flow'] += pi[1]

    data[str(i)] = {'ttt' : ttt, 'Pi_ue' : Pi_ue, 'Pi_so' : Pi_so, 'links': output_links}

    elog.append((ttt))

with open('../network/data.pickle', 'wb') as fw:
    pickle.dump(data, fw)

plt.plot(elog)

plt.legend()
plt.show()


# print(optimal_flow)
