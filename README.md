### Static Mixed Traffic Assignment
## [CEE301] 2024년 교통계획 및 실험 텀프로젝트
"CAV 혼재 상황(Mixed Traffic)에서의 교통혼잡세 추정 방안 연구"에서 사용한 Static Mixed Traffix Assignment 모델입니다.
assignment파일에 있는 bilevel.py파일을 실행하면 data.pickle에 데이터가 저장됩니다.

e는 0부터 1까지 10%씩 증가하도록 되어있습니다.

코드 정리를 안해서 변수명도 마음대로고 보기 조금 더럽습니다.
주석은 보시는 분들 이해하기 쉬우라고 달아놓은 것이 아니라,
저도 제 코드 보면서 어지러워서 써놓은 거라 매우 불친절합니다.

# data.pickle 데이터 형식
data.pickle은 



    with open('data.pickle', 'rb') as fr:
    
        data = pickle.load(fr)
    
로 다시 읽을 수 있습니다.



data 형식

data[e(string)] => e별로 데이터가 나눠져 있습니다.
Ex)
data['1'] -> CAV 비율이 10%일 때 데이터
data['10'] -> CAV 비율이 100%일 때 데이터

data[e]에는 4가지 데이터가 들어있습니다.

1. data[e]['ttt'] : Total Travel Time(float)

2. data[e]['Pi_ue'] : HDV차량들의 경로 집합

od별로 HDV차량들의 경로 집합이 저장됩니다.

Ex)
data[e]['Pi_ue'][('2','17')] = [[[path1], flow1], [[path2], flow2], ...]

=> 2에서 17로 가는 HDV 차량들의 경로(path)와 각 경로로 이동하는 차량 수(flow)

3. data[e]['Pi_so'] : CAV차량들의 경로 집합
   
od별로 CAV차량들의 경로 집합이 저장됩니다.

Ex)
data[e]['Pi_ue'][('2','17')] = [[[path1], flow2], [[path2], flow2], ...]

=> 2에서 17로 가는 CAV 차량들의 경로(path)와 각 경로로 이동하는 차량 수(flow)

4. data[e]['links'] : e일때 각 도로들 정보

data[e]['links'][link(tuple)]['capacity'] : link에서의 도로 용량

data[e]['links'][link(tuple)]['flow'] : link에서의 통행량

data[e]['links'][link(tuple)]['travel_times'] : BPR function을 통해 구한 travel time

data[e]['links'][link(tuple)]['marginal_travel_times'] : link의 marginal travel time


Ex)
data[e]['links'][('2','17')]['travel_times'] -> 2에서 17로 연결된 도로에서의 travel time
