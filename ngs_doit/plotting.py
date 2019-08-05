# type: ignore
import csv
import plotly.graph_objects as go
import itertools
#_, poss, depths = itertools.tee(csv.DictReader(open('depth.hist'),  fieldnames=('chrom', 'pos', 'depth')))
print(list(csv.reader(open('depth.hist'), delimiter='\t'))[:4])
#_, poss, depths = itertools.tee(csv.reader(open('depth.hist'), delimiter='\t'))
reader = csv.reader(open('depth.hist'), delimiter='\t')
_, depthX, depthY = zip(*list(reader))
rows = csv.DictReader(open('qual.hist'), delimiter='\t')
def avg(*args): return sum(map(float, args)) / len(args)
#qualX, qualY = (x[0], avg(x[1] + x[3])
qualX, qualY = zip(*[(x['#BaseNum'], avg(x['Read1_linear'], x['Read2_linear'] ) ) for x in rows])
#_,  xs, ys = itertools.tee(reader, 3)
#xs, ys = (int(x[1]) for x in xs), (int(y[2]) for y in ys)
from plotly.subplots import make_subplots

fig = make_subplots(specs=[[{"secondary_y": True}]])
#fig = go.Figure()
fig.add_trace(go.Scatter(x=list(depthX), y=list(depthY), fill='tozeroy', stackgroup='one', 
    mode='lines',
    line=dict(width=0.5, color='rgb(131, 90, 241)'))) # fill down to xaxis
fig.add_trace(go.Scatter(x=qualX, y=qualY, fill='tonexty', stackgroup='one',
    mode='lines',
    line=dict(width=0.5, color='rgb(111, 231, 219)')), secondary_y=True) # fill down to xaxis
#fig.update_layout(yaxis_range=(0, 100))

fig.show()
