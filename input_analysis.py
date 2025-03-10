from operate_analysis import analyze

'''

This script can be run to get custom feedback from your openAI account on a 
reaction mechanism. You can draw your reaction mechanism in ChemDoodle. To run
the analysis you have to copy your openai key, the exercise number and the
ChemDoodle drawing as JSON code. 

'''

# Copy the openAI key between quotation marks
# find your key under "specify later"
ai_key = False

# integer number which corresponds to the exercise number of the reaction mechanism
exercise = 3


# copy json_code from your reaction mechanism in ChemDoodle
json_code = '''{"m":[{"a":[{"x":190.3717291886041,"y":265.0121472612378,"i":"a0"},{"x":207.69223726429283,"y":255.01214726123783,"i":"a1"},{"x":225.01274533998162,"y":265.0121472612378,"i":"a2"},{"x":207.69223726429283,"y":235.01214726123777,"i":"a3"},{"x":225.01274533998156,"y":245.01214726123777,"i":"a4","l":"Cl"},{"x":190.37172918860404,"y":245.01214726123783,"i":"a5"},{"x":207.69223726429283,"y":275.0121472612378,"i":"a6"}],"b":[{"b":0,"e":1,"i":"b0"},{"b":1,"e":2,"i":"b1"},{"b":1,"e":3,"i":"b2"},{"b":1,"e":4,"i":"b3"},{"b":1,"e":5,"i":"b4"},{"b":1,"e":6,"i":"b5"}]},{"a":[{"x":461.27110749852955,"y":247.81870698697776,"i":"a7","l":"Cl","c":-1}]},{"a":[{"x":370.80377547566275,"y":269.5086908690869,"i":"a8"},{"x":388.1242835513515,"y":259.5086908690869,"i":"a9"},{"x":405.44479162704033,"y":269.5086908690869,"i":"a10"},{"x":388.1242835513515,"y":239.50869086908682,"i":"a11"},{"x":405.44479162704033,"y":249.50869086908688,"i":"a12","l":"O","c":1},{"x":422.76529970272895,"y":259.5086908690869,"i":"a13","l":"H"},{"x":370.8037754756627,"y":249.50869086908688,"i":"a14"}],"b":[{"b":0,"e":1,"i":"b6"},{"b":1,"e":2,"i":"b7"},{"b":1,"e":3,"i":"b8"},{"b":1,"e":4,"i":"b9"},{"b":4,"e":5,"i":"b10"},{"b":1,"e":6,"i":"b11"}]},{"a":[{"x":614.2691220103162,"y":260.9938393839383,"i":"a15"},{"x":631.5896300860051,"y":250.9938393839384,"i":"a16"},{"x":648.9101381616938,"y":260.9938393839383,"i":"a17"},{"x":631.5896300860051,"y":230.99383938393837,"i":"a18"},{"x":648.9101381616938,"y":240.99383938393837,"i":"a19","l":"O"}],"b":[{"b":0,"e":1,"i":"b12"},{"b":1,"e":2,"i":"b13"},{"b":1,"e":3,"i":"b14"},{"b":1,"e":4,"i":"b15"}]},{"a":[{"x":139.4187802953374,"y":252.17645764576457,"i":"a20","l":"O"}]},{"a":[{"x":549.834621879496,"y":244.50869086908688,"i":"a21","l":"Cl"}]}],"s":[{"i":"s0","t":"Line","x1":275.2554501188366,"y1":254.54703098216805,"x2":322.9440469895573,"y2":254.54703098216805,"a":"synthetic"},{"i":"s1","t":"Pusher","o1":"b3","o2":"a4","e":2},{"i":"s2","t":"Pusher","o1":"a7","o2":"a13","e":2},{"i":"s3","t":"Pusher","o1":"b10","o2":"a12","e":2},{"i":"s4","t":"Line","x1":484.6813696781835,"y1":248.96479647964793,"x2":526.9035919004058,"y2":248.96479647964793,"a":"synthetic"},{"i":"s5","t":"Pusher","o1":"a20","o2":"a1","e":2}]}'''

''' now run the script '''
output = analyze(ai_key,exercise,json_code)

print(output)