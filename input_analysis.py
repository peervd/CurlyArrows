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
exercise =  4  # paste integer number between 1 and 9

# copy json_code from your reaction mechanism in ChemDoodle
json_code = '''
{"m":[{"a":[{"x":39.22304240956336,"y":225.9619625404376,"i":"a0"},{"x":56.543550485252126,"y":215.9619625404376,"i":"a1"},{"x":73.86405856094092,"y":225.9619625404376,"i":"a2"},{"x":56.5435504852521,"y":195.9619625404376,"i":"a3"},{"x":91.18456663662971,"y":215.9619625404376,"i":"a4"},{"x":108.50507471231845,"y":225.96196254043764,"i":"a5"},{"x":91.18456663662971,"y":195.9619625404376,"i":"a6"}],"b":[{"b":0,"e":1,"i":"b0"},{"b":1,"e":2,"i":"b1"},{"b":1,"e":3,"i":"b2"},{"b":2,"e":4,"i":"b3"},{"b":4,"e":5,"i":"b4","o":2},{"b":4,"e":6,"i":"b5"}]},{"a":[{"x":152.51899184195412,"y":213.5850941150504,"i":"a7","l":"H","c":1}]},{"a":[{"x":251.72100016672982,"y":228.45782893653188,"i":"a8"},{"x":269.0415082424187,"y":218.45782893653185,"i":"a9"},{"x":286.3620163181075,"y":228.45782893653188,"i":"a10"},{"x":269.0415082424187,"y":198.45782893653185,"i":"a11"},{"x":303.68252439379614,"y":218.45782893653185,"i":"a12","c":1},{"x":303.68252439379614,"y":198.45782893653185,"i":"a13"},{"x":321.0030324694852,"y":228.45782893653188,"i":"a14"}],"b":[{"b":0,"e":1,"i":"b6"},{"b":1,"e":2,"i":"b7"},{"b":1,"e":3,"i":"b8"},{"b":2,"e":4,"i":"b9"},{"b":4,"e":5,"i":"b10"},{"b":4,"e":6,"i":"b11"}]},{"a":[{"x":363.87782187548976,"y":220.24956172872044,"i":"a15","l":"O"},{"x":381.19832995117844,"y":210.24956172872044,"i":"a16"}],"b":[{"b":0,"e":1,"i":"b12"}]},{"a":[{"x":486.88150024908754,"y":229.40989573848472,"i":"a17"},{"x":504.2020083247764,"y":219.40989573848472,"i":"a18"},{"x":521.5225164004651,"y":229.40989573848472,"i":"a19"},{"x":504.2020083247764,"y":199.40989573848475,"i":"a20"},{"x":538.8430244761539,"y":219.40989573848472,"i":"a21"},{"x":538.8430244761539,"y":199.40989573848475,"i":"a22"},{"x":556.1635325518428,"y":229.40989573848478,"i":"a23"},{"x":538.8430244761539,"y":239.40989573848472,"i":"a24","l":"O","c":1},{"x":521.5225164004653,"y":249.40989573848472,"i":"a25"},{"x":556.1635325518428,"y":249.40989573848472,"i":"a26","l":"H"}],"b":[{"b":0,"e":1,"i":"b13"},{"b":1,"e":2,"i":"b14"},{"b":1,"e":3,"i":"b15"},{"b":2,"e":4,"i":"b16"},{"b":4,"e":5,"i":"b17"},{"b":4,"e":6,"i":"b18"},{"b":4,"e":7,"i":"b19"},{"b":7,"e":8,"i":"b20"},{"b":7,"e":9,"i":"b21"}]},{"a":[{"x":587.1968553513519,"y":422.1865052943491,"i":"a27","l":"S"},{"x":604.5173634270406,"y":412.186505294349,"i":"a28","l":"O"},{"x":569.8763472756633,"y":412.186505294349,"i":"a29","l":"O","c":-1},{"x":587.1968553513519,"y":442.186505294349,"i":"a30","l":"O"},{"x":569.8763472756633,"y":432.186505294349,"i":"a31","l":"O"}],"b":[{"b":0,"e":1,"i":"b22","o":2},{"b":0,"e":2,"i":"b23"},{"b":0,"e":3,"i":"b24"},{"b":0,"e":4,"i":"b25","o":2}]},{"a":[{"x":765.623371442933,"y":228.32067990231394,"i":"a32"},{"x":782.9438795186219,"y":218.3206799023139,"i":"a33"},{"x":800.2643875943105,"y":228.32067990231394,"i":"a34"},{"x":782.9438795186219,"y":198.32067990231397,"i":"a35"},{"x":817.5848956699994,"y":218.3206799023139,"i":"a36"},{"x":817.5848956699994,"y":198.32067990231397,"i":"a37"},{"x":834.9054037456882,"y":228.32067990231394,"i":"a38"},{"x":817.5848956699994,"y":238.32067990231397,"i":"a39","l":"O"},{"x":800.2643875943107,"y":248.32067990231397,"i":"a40"}],"b":[{"b":0,"e":1,"i":"b26"},{"b":1,"e":2,"i":"b27"},{"b":1,"e":3,"i":"b28"},{"b":2,"e":4,"i":"b29"},{"b":4,"e":5,"i":"b30"},{"b":4,"e":6,"i":"b31"},{"b":4,"e":7,"i":"b32"},{"b":7,"e":8,"i":"b33"}]}],"s":[{"i":"s0","t":"Line","x1":179.17686229663428,"y1":216.441294520909,"x2":213.78026688488217,"y2":216.441294520909,"a":"synthetic"},{"i":"s1","t":"Pusher","o1":"b4","o2":"a7","e":2},{"i":"s2","t":"Pusher","o1":"a15","o2":"a12","e":2},{"i":"s3","t":"Line","x1":417.19356278485026,"y1":221.2016285306733,"x2":448.6117672492947,"y2":220.24956172872044,"a":"synthetic"},{"i":"s4","t":"Line","x1":693.1601709541086,"y1":218.68761107106815,"x2":737.9073106458931,"y2":218.68761107106815,"a":"synthetic"},{"i":"s5","t":"Pusher","o1":"a29","o2":"a26","e":2},{"i":"s6","t":"Pusher","o1":"b21","o2":"a24","e":2}]}

'''


''' now run the script '''
output = analyze(ai_key,exercise,json_code)
print(output)