import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt


# Define the lists
list1 = np.array([(1, 0.5), (2, 0.3), (3, 0.2)])
list2 = np.array([(7, 0.1), (2, 0.2), (3, 0.3), (4, 0.4)])
list5 = np.array([(1, 0.6), (2, 0.3), (3, 0.1), (4, 0.4), (5.5, 0.5)])

# Create a dictionary to store the lists
lists = {
    "list1": list1,
    "list2": list2,
    "list5": list5,
    # Add more lists here...
}

# Determine the common X-axis range
x_values = np.unique(np.concatenate([lst[:, 0] for lst in lists.values()]))

# Pad the lists based on the common X-axis range
padded_lists = {}

# Pad and fill the lists with their respective X-axis values
for name, lst in lists.items():
    padded_lst = np.zeros((len(x_values), 2))
    for i, x in enumerate(x_values):
        if x in lst[:, 0]:
            index = np.where(lst[:, 0] == x)[0][0]
            padded_lst[i] = lst[index]
    for i, x in enumerate(x_values):
        if padded_lst[i, 0] == 0:
            padded_lst[i, 0] = x
    padded_lists[name] = padded_lst

# Print the padded lists
for name, padded_lst in padded_lists.items():
    print(f"Padded {name}:")
    print(padded_lst)
    print()


# ----------------------------------------------
# def _cal_cosine(list1, list2):

# Define the cosine similarity score function
def cosine_score(weights, lists):
    weights = np.array(weights, dtype=float)
    print("weights, ", weights)
    combined_list = np.sum([w * np.array(lst) for w, lst in zip(weights, lists[:-1])], axis=0)
    print("prediction_list: ", lists)
    print("combined_list: ", combined_list)
    dot_product = np.dot(combined_list[:, 1], lists[-1][:, 1])
    print("dot_product: ", dot_product)
    magnitude = np.linalg.norm(combined_list[:, 1]) * np.linalg.norm(lists[-1][:, 1])
    cosine_score = dot_product / magnitude
    return -cosine_score  # Negate the score for minimization

# Define the lists of values
list1 = np.array([(15.023475, 0.1263170347145845), (16.026829, 0.01308932482555683), (16.029751, 0.0023420407534433746), (26.01565, 0.012069344065364484), (27.019004, 0.0017713988703897193), (27.021926, 0.0001809640881740451), (28.022358, 0.00012200596134020557), (29.002739, 0.08149537867007935), (30.006093, 0.008815630910652338), (30.006956, 0.0013232355589213839), (30.009015, 0.0010380313602479376), (31.006984, 0.0035009522148635387), (32.010338, 0.00036699950695780637), (39.023475, 0.2078256687993354), (40.026829, 0.03714962024949196), (40.029751, 0.0037780631069858563), (40.0313, 0.008274437489666109), (41.030183, 0.004328315079511918), (41.033105, 0.0009445157767464649), (41.034654, 0.0014963896621793286), (41.037576, 0.0001724054938468704), (42.038008, 0.00018042859947706472), (42.04093, 5.3205490044012106e-05), (42.04695, 0.0054464119182829175), (43.018389, 0.17963458679475755), (43.050304, 0.0009792046866610247), (43.053226, 0.00014006515118522665), (44.021743, 0.02616641695156233), (44.022606, 0.003898776218657787), (44.024665, 0.003251804077940207), (44.053658, 0.0001050488655120529), (44.05658, 2.4760254552775607e-05), (45.022634, 0.007495038210097604), (45.025097, 0.002438853107746798), (45.02596, 0.0008129510194850524), (45.028019, 0.0008129510194850524), (46.025988, 0.0014080724698130257), (49.007825, 0.1301986463289096), (50.011179, 0.027213730013324076), (50.014101, 0.0015745631264442285), (51.014533, 0.003716578884017371), (51.017455, 0.0005951289068414936), (52.0313, 0.17419512379150132), (53.002739, 0.2233202635724143), (53.034654, 0.036533435347248865), (53.037576, 0.0035619074961799452), (54.006093, 0.04067692608001754), (54.006956, 0.00406642206947142), (54.009015, 0.0026896854503044884), (54.010564, 0.04550588926632016), (54.038008, 0.004973933198632075), (54.04093, 0.0013795208413436284), (55.006984, 0.009908645116278307), (55.009447, 0.004191571972542914), (55.013918, 0.008132481020640856), (55.014781, 0.0009928965981392487), (55.01684, 0.0006546966133808274), (56.010338, 0.0026896854503044884), (56.01326, 0.0010166055173678548), (56.014809, 0.0019087522659332191), (56.017272, 0.0009710720227902142), (56.018135, 0.00020703324746820592), (56.020194, 0.00020703324746820592), (56.026214, 0.028360398844580848), (57.018163, 0.0003585921034709183), (57.029568, 0.005127138086085856), (57.030431, 0.0005627765363195793), (57.03249, 0.0005773965267420138), (58.030459, 0.0013039449406510794), (58.032922, 0.0006191888077562544), (58.033785, 0.000129109788376746), (58.035844, 0.00018258881375750444), (59.033813, 0.000258219576753485), (59.036735, 0.000129109788376746), (62.01565, 0.12392101217560454), (63.019004, 0.028839172286557684), (63.021926, 0.0018006906090031782), (64.022358, 0.004555426948580595), (64.02528, 0.0005694283685725741), (65.002739, 0.011292387382645834), (66.006093, 0.002359526336711806), (66.006956, 0.00024237238268959234), (66.009015, 0.00013671646356364203), (66.010564, 0.10845677736920131), (67.006984, 0.00047359971404225236), (67.009447, 0.00032270381100924775), (67.01031, 7.307802178480551e-05), (67.012369, 5.1673964759734223e-05), (67.013918, 0.02243020045357853), (67.014781, 0.002045583338431583), (67.01684, 0.001568890785559512), (67.018389, 0.17142664078801853), (68.010338, 0.00010334792951948596), (68.014809, 0.00488628262882539), (68.017272, 0.002976761029614616), (68.020194, 0.0004961268282418879), (68.021743, 0.03566810715380987), (68.022606, 0.0034196891224795), (68.024665, 0.0031381218684416088), (68.026214, 0.13063065214219263), (68.997653, 0.13052941117830602), (69.018163, 0.001215257576898053), (69.022634, 0.007884433743236755), (69.025097, 0.005144510286756046), (69.02596, 0.0007845304671104025), (69.028019, 0.0007845304671104025), (69.029568, 0.027384901757101916), (69.030431, 0.0019829181339871137), (69.03249, 0.0026737662222058468), (69.034039, 0.15144284160158905), (70.001007, 0.023787702901988738), (70.00187, 0.0037150150286066955), (70.003929, 0.001573900584590126), (70.025988, 0.0019216993320878256), (70.030459, 0.005479590215189235), (70.032922, 0.0037337113870233034), (70.035844, 0.0010355452050255167), (70.037393, 0.031398328296627726), (70.038256, 0.002857172804762116), (70.040315, 0.0033233491516643943), (71.001898, 0.008025349697029418), (71.004361, 0.0024527468731297343), (71.005224, 0.0005948784895652702), (71.033813, 0.0014644820733975131), (71.038284, 0.006998615378777513), (71.040747, 0.0040996432139259065), (71.04161, 0.000980002176562751), (71.043669, 0.000980002176562751), (72.005252, 0.0016825704558144493), (72.008174, 0.0005948784895652702), (72.021128, 0.035012267596300994), (72.041638, 0.0015495194949303875), (73.007825, 0.12218010134927049), (73.024482, 0.00632595094002378), (73.025345, 0.000873896861300095), (73.027404, 0.0007135337733767904), (74.011179, 0.031401201934608186), (74.014101, 0.0014940858396303314), (74.025373, 0.0021818271885646265), (74.027836, 0.0007651797438476834), (74.028699, 0.0001595510021512416), (74.030758, 0.0002256391911325047), (75.014533, 0.00529746214105947), (75.017455, 0.0005647113523320837), (75.028727, 0.00047865301612778455), (75.031649, 0.0001595510021512416), (76.017887, 0.0005647113523320837), (76.0313, 0.12881082568229038), (77.002739, 0.021811219214385804), (77.034654, 0.03257111969716383), (77.037576, 0.0026603076656354615), (77.039125, 0.17906123073816402), (78.006093, 0.005091360673835247), (78.006956, 0.00041376529084715026), (78.009015, 0.00026550861822788), (78.010564, 0.11059225749247528), (78.038008, 0.005611925317663016), (78.04093, 0.0010303327284794065), (78.042479, 0.045372806171286535), (78.045401, 0.0039662156238958884), (79.006984, 0.0009832528361224644), (79.009447, 0.000770824671393895), (79.012369, 0.0001003528223601779), (79.013918, 0.0257241489101192), (79.014781, 0.002217876656427309), (79.01684, 0.0016090167936990277), (79.018389, 0.1299027364930153), (79.041362, 0.0005948628781424679), (79.045833, 0.007579700394757615), (79.048755, 0.0011695731524894182), (80.010338, 0.00024581320903061286), (80.014809, 0.005062653140955789), (80.017272, 0.00407052628924015), (80.018135, 0.0005088157861550182), (80.020194, 0.0005088157861550182), (80.021743, 0.03030372411130192), (80.022606, 0.0019821718845780176), (80.024665, 0.0023905892008820083), (80.049187, 0.0011695731524894182), (80.997653, 0.07006111077059549), (81.013692, 0.0001003528223601779), (81.018163, 0.0014391483710599114), (81.022634, 0.005378825762379749), (81.025097, 0.004855309753232855), (81.028019, 0.0005976473002205019), (81.034039, 0.12751429753148572), (82.001007, 0.014635557570421592), (82.00187, 0.002004538199932825), (82.003929, 0.0008492412064052307), (82.025988, 0.0017929419368986435), (82.028451, 0.0005976473002205019), (82.037393, 0.02945914277272932), (82.038256, 0.0021146240074255917), (82.040315, 0.002812712897051738), (83.001898, 0.00431838467154351), (83.004361, 0.002004538199932825), (83.005224, 0.0004539385072182879), (83.007283, 0.0003209829966957499), (83.038284, 0.005894163384054195), (83.040747, 0.0045429412158603895), (83.04161, 0.0005864911677244617), (83.043669, 0.0008294237636079652), (84.005252, 0.0010150373597448672), (84.041638, 0.0013114341192351257), (84.044101, 0.0005864911677244617), (87.023475, 0.044002065645769715), (88.026829, 0.01200136779080091), (88.029751, 0.000817062608448204), (89.030183, 0.002246922217143852), (89.033105, 0.00020426565211205318), (90.033537, 0.0005003466196513876), (92.026214, 0.024859309248990404), (92.997653, 0.028864150383175434), (93.029568, 0.006287282170732986), (93.030431, 0.00041443252291240994), (93.03249, 0.0005140402726532032), (93.034039, 0.12919737075071258), (94.001007, 0.006735655399708457), (94.00187, 0.0007866572370615656), (94.003929, 0.00035180382238255), (94.005478, 0.08239185387948997), (94.030459, 0.0011551618600440456), (94.032922, 0.0010843691719062249), (94.033785, 0.00011494289928251044), (94.035844, 0.00019908694152658034), (94.037393, 0.03273956204456818), (94.038256, 0.003104912689733811), (94.040315, 0.002865704736277812), (95.001898, 0.0018568158206250479), (95.004361, 0.0010213569246569458), (95.005224, 0.0001329693428741641), (95.007283, 0.0001329693428741641), (95.008832, 0.019149520003810743), (95.009695, 0.00207828000724157), (95.011754, 0.0011998954779803943), (95.033813, 0.00025702013632659807), (95.036276, 0.00011494289928251044), (95.038284, 0.006064376275629991), (95.040747, 0.005476551297402148), (95.04161, 0.00084505020410869), (95.043669, 0.00084505020410869), (96.005252, 0.00047942779050583005), (96.006115, 0.0001329693428741641), (96.009723, 0.005090725395747787), (96.012186, 0.0030355221316435686), (96.013049, 0.00037944026645544705), (96.015108, 0.00037944026645544705), (96.021128, 0.04880869962405777), (96.041638, 0.0014636698884627013), (96.044101, 0.00084505020410869), (97.008606, 0.0001329693428741641), (97.013077, 0.0015644723139480695), (97.024482, 0.01129314448142212), (97.025345, 0.001251558215330705), (97.027404, 0.0010052757321043439), (98.025373, 0.0031948163071164213), (98.027836, 0.0017411887069647984), (98.028699, 0.0003178960989942985), (98.030758, 0.00038934161687865766), (99.028727, 0.0007455324486829468), (101.039125, 0.12787948900249513), (102.010564, 0.13125432640269297), (102.042479, 0.037810509677658054), (102.045401, 0.002865333612912988), (103.013918, 0.03589884488951018), (103.014781, 0.002440489480636429), (103.01684, 0.0019293763411231276), (103.018389, 0.20757573521761555), (103.045833, 0.007366023034372067), (103.048755, 0.0008449407658015173), (104.014809, 0.00603990342621094), (104.017272, 0.006542833658800891), (104.018135, 0.0010567629439499195), (104.020194, 0.0006101223701591085), (104.021743, 0.05663574652628242), (104.022606, 0.003478601382446144), (104.024665, 0.0038591616848661944), (104.026214, 0.2622263255649338), (104.049187, 0.001463480335754385), (105.018163, 0.0017256866611724546), (105.020626, 0.001494488487547294), (105.022634, 0.009696023890836217), (105.025097, 0.010612694840784458), (105.02596, 0.000964790421216549), (105.028019, 0.000964790421216549), (105.029568, 0.0717271644235165), (105.030431, 0.006335264322391574), (105.03249, 0.005452527417721352), (105.034039, 0.050667680414009306), (106.021517, 0.0006101223701591085), (106.025988, 0.002157336965880858), (106.028451, 0.002363244240705398), (106.030459, 0.01237375427128366), (106.032922, 0.013244143365764736), (106.033785, 0.0017242405644515807), (106.035844, 0.0021117547883574127), (106.037393, 0.013963915480443357), (106.038256, 0.0010273784440437425), (106.040315, 0.0011303620403972908), (107.033813, 0.002986472262101378), (107.036276, 0.0017242405644515807), (107.038284, 0.002415173215394958), (107.040747, 0.002549450373075963), (107.04161, 0.0004082387859011444), (107.043669, 0.00033332557289036854), (108.041638, 0.0005270340063570214), (108.044101, 0.0004713935458673409), (108.057514, 0.1279606142079484), (109.060868, 0.03539062712757715), (109.061731, 0.002595333899759175), (109.06379, 0.003812486938646202), (110.061759, 0.005864112557671752), (110.064222, 0.0062162673986937145), (110.065085, 0.0010312810886663431), (110.067144, 0.0005954104141516825), (111.065113, 0.0015753079247985903), (111.067576, 0.000842037482871493), (111.994914, 0.12775284045073593), (112.068467, 0.0005954104141516825), (112.998268, 0.03794492719240261), (112.999131, 0.001981936784086308), (113.999159, 0.0053113754880814105), (114.001622, 0.007535151975456918), (114.04695, 0.127179794392449), (115.002513, 0.0018897025466273875), (115.004976, 0.0011951528295126165), (115.050304, 0.039997814081483216), (115.053226, 0.003327565221532451), (116.005867, 0.0005975764147563068), (116.053658, 0.008388394590840271), (116.05658, 0.0008452023965938745), (117.057012, 0.0015812289363563783), (118.041864, 0.07164191288115532), (119.045218, 0.021245847345865), (119.046081, 0.0013821259021092576), (119.04814, 0.0018962610733145346), (120.021128, 0.012209191922157168), (120.046109, 0.003215267159628882), (120.048572, 0.004036520333836018), (120.049435, 0.0005806090060903714), (120.051494, 0.0003352147659601988), (121.024482, 0.003339660675839074), (121.025345, 0.00038551907262434547), (121.027404, 0.000254203742578921), (121.049463, 0.0011117816212614917), (121.051926, 0.0008211051308489662), (122.025373, 0.0008158316829467369), (122.027836, 0.0006148368212848461), (122.028699, 0.0001136833697082682), (122.030758, 9.845268615517859e-05), (123.028727, 0.00020494561219298418), (123.03119, 8.03862816288538e-05), (123.032053, 5.68416848541341e-05), (130.041864, 0.022243494157291287), (131.045218, 0.00699497140468396), (131.046081, 0.00044403194917840255), (131.04814, 0.0005827185091709393), (132.021128, 0.18053495010550358), (132.046109, 0.0010254478408224638), (132.048572, 0.0014652306751110247), (132.049435, 0.00020931866398781338), (132.051494, 0.00014801064673469123), (133.024482, 0.053126764000715664), (133.025345, 0.004781177086934385), (133.027404, 0.0037798523728303835), (133.028953, 0.18048974818140776), (133.049463, 0.00029602129346938247), (133.051926, 0.000276902572033497), (134.025373, 0.01171144423287573), (134.027836, 0.010454555259574137), (134.028699, 0.001463930529109466), (134.030758, 0.001463930529109466), (134.032307, 0.05333486814547309), (134.03317, 0.00462936296232387), (134.035229, 0.004053447177300841), (134.052817, 0.00010465933199391534), (135.028727, 0.003273448296513353), (135.03119, 0.001690401370112523), (135.032053, 0.001195294271733596), (135.033198, 0.011308034157245983), (135.035661, 0.010386017852399947), (135.036524, 0.0014639330543991939), (135.038583, 0.0011952963336240244), (136.036552, 0.0037798588931004223), (136.037415, 0.0008452021430329655), (136.039015, 0.00207031397993761), (136.039878, 0.0008452021430329655), (144.057514, 0.21620813015519058), (145.060868, 0.07149486357778058), (145.061731, 0.004090014205007777), (145.06379, 0.006466880275133713), (145.065339, 0.37889314120879647), (146.061759, 0.01012226891616668), (146.064222, 0.015234959662963133), (146.065085, 0.0020450071025038884), (146.067144, 0.0014460383897551522), (146.068693, 0.12609392810569378), (146.069556, 0.008018296376698028), (146.071615, 0.0116196154663811), (147.044603, 0.12692118876655423), (147.065113, 0.00354205620337595), (147.067576, 0.002504611960751011), (147.069584, 0.016432617782116327), (147.072047, 0.027013449571627605), (147.074969, 0.004009148188349014), (148.047957, 0.03983675996023798), (148.04882, 0.0034332294893998014), (148.050879, 0.0035858912248818335), (148.072938, 0.0050712159009282465), (148.075401, 0.00537883683732275), (148.07586, 0.0017929455762037054), (149.048848, 0.007883525544116622), (149.051311, 0.008150839196299957), (149.052174, 0.0008452026501547073), (149.054233, 0.001035157611063558), (150.052202, 0.0024641680235903754), (150.054665, 0.0011952970508024693), (151.055556, 0.0008452026501547073)])




list2 = np.array([(51.026928, 0.03012272320350666), (28.026829, 0.004391425701048351), (28.029751, 0.0005446893136763442), (29.030183, 0.000408516993513792), (29.033105, 0.00013617232841908134), (39.023475, 0.12639400529867006), (40.026829, 0.022593404009162908), (40.029751, 0.002297716788892669), (41.002739, 0.006660620772566236), (41.030183, 0.0026323653004690453), (41.033105, 0.0005744291972231664), (41.039125, 0.12294526171155472), (42.006093, 0.0010013217706568264), (42.006956, 0.0001243676188644567), (42.009015, 8.531548177831184e-05), (42.042479, 0.02245800472027476), (42.045401, 0.002681747615646734), (43.006984, 0.0003134690941630608), (43.009447, 6.032715570565869e-05), (43.018389, 0.12104174945898677), (43.045833, 0.0025624985361980177), (43.048755, 0.0007908042099667905), (44.010338, 3.0163577852812294e-05), (44.021743, 0.01763150928450725), (44.022606, 0.0026270814695313676), (44.024665, 0.0021911373612113153), (45.022634, 0.005050322175699423), (45.025097, 0.0016433530541222837), (45.02596, 0.0005477843403028281), (45.028019, 0.0005477843403028281), (46.025988, 0.0009487903089950982), (49.007825, 0.014075155354570011), (50.011179, 0.0029419466984873747), (50.014101, 0.00017021851797363022), (51.014533, 0.00040178163640740776), (51.017455, 6.433655077048819e-05), (52.0313, 0.11892032598366069), (53.034654, 0.024940813188306887), (53.037576, 0.0024316593447033658), (54.010564, 0.004550181543821028), (54.038008, 0.003395627526923557), (54.04093, 0.0009417776145685771), (55.013918, 0.0008131752975758783), (55.014781, 9.928077109613443e-05), (55.01684, 6.546380029126692e-05), (55.018389, 0.17205821294702817), (56.014809, 0.00019085813885201625), (56.017272, 9.709850894159682e-05), (56.018135, 2.0701471321040018e-05), (56.020194, 2.0701471321040018e-05), (56.021743, 0.030742842895019808), (56.022606, 0.003227946274182001), (56.024665, 0.0031315678352621537), (56.026214, 0.01982717966353146), (57.018163, 3.585600011952868e-05), (57.022634, 0.007750235126061427), (57.025097, 0.003587661693320259), (57.028019, 0.0007828919588155396), (57.029568, 0.0035844590391571245), (57.030431, 0.0003934455067849072), (57.03249, 0.0004036665610928271), (58.025988, 0.00175059963894951), (58.030459, 0.0009116074407600209), (58.032922, 0.0004328841707872145), (58.033785, 9.026258708470763e-05), (58.035844, 0.0001276505748300922), (59.033813, 0.00018052517416942667), (59.036735, 9.026258708470763e-05), (63.023475, 0.1253191150062254), (64.026829, 0.029257108382929047), (64.029751, 0.0023039949669369647), (65.030183, 0.004679436027838305), (65.033105, 0.000575998741734241), (66.010564, 0.10938463562190365), (66.033537, 0.000575998741734241), (67.013918, 0.022622092994602255), (67.014781, 0.0020630835023511974), (67.01684, 0.0015823127984413308), (67.018389, 0.17966043378777777), (68.014809, 0.004928085250774574), (68.017272, 0.0030022275090237204), (68.020194, 0.0005003712413909537), (68.021743, 0.03738128201185813), (68.022606, 0.003583940211041883), (68.024665, 0.003288849000198903), (68.026214, 0.12678493675213384), (69.018163, 0.0012256542233708267), (69.022634, 0.008263130981097265), (69.025097, 0.005391606260821484), (69.02596, 0.0008222122500497262), (69.028019, 0.0008222122500497262), (69.029568, 0.026578700942702704), (69.030431, 0.0019245417984177224), (69.03249, 0.0025950515886835114), (70.025988, 0.0020140004728874806), (70.030459, 0.005318273218939144), (70.032922, 0.0036237923817389723), (70.035844, 0.0010050591585520337), (71.033813, 0.0014213682930115767), (73.007825, 0.12461385884222513), (74.011179, 0.03202669585425775), (74.014101, 0.0015238471720172065), (74.01565, 0.17667116831745544), (75.014533, 0.005402984546402639), (75.017455, 0.0005759600783514371), (75.019004, 0.04515818658504814), (75.021926, 0.0025815216624137734), (75.023475, 0.12485659746023062), (76.017887, 0.0005759600783514371), (76.022358, 0.007787471601750109), (76.02528, 0.0008163488282291822), (76.026829, 0.031581445495165866), (76.029751, 0.002306382676153405), (77.025712, 0.001154491584509092), (77.030183, 0.005439592807060664), (77.033105, 0.0005765956690383513), (77.039125, 0.12498822434943214), (78.033537, 0.0011531913380767025), (78.042479, 0.031671101855614385), (78.045401, 0.002768495704927965), (79.045833, 0.0052907784087933615), (79.048755, 0.0008163848253125147), (80.049187, 0.0008163848253125147), (86.01565, 0.12328471529049746), (87.019004, 0.03373111619640911), (87.021926, 0.0018099989566941392), (88.022358, 0.006164634349760552), (88.02528, 0.0005723719265681947), (88.0313, 0.1243550424582294), (89.025712, 0.0014020191631858407), (89.034654, 0.034006561809530246), (89.037576, 0.0025821149194915237), (90.010564, 0.12477431531172023), (90.038008, 0.006271935482521122), (90.04093, 0.0010000488081184369), (91.013918, 0.031871869516649706), (91.014781, 0.0019140530075622169), (91.01684, 0.0018249779063652621), (91.041362, 0.0008165364325895514), (92.014809, 0.005161817012426019), (92.017272, 0.005505266071001078), (92.020194, 0.0005771086863599715), (93.018163, 0.0018249779063652621), (93.020626, 0.0008161549312135949), (101.039125, 0.12358148486430151), (102.042479, 0.03653970598326251), (102.045401, 0.0027690303212618526), (103.018389, 0.14168777073038524), (103.045833, 0.007118452468281182), (103.048755, 0.0008165424750648966), (104.021743, 0.038658625780838), (104.022606, 0.0023744358878060394), (104.024665, 0.0026342000689221485), (104.026214, 0.12410211155727713), (104.049187, 0.001414293053350444), (105.022634, 0.0066183458707294715), (105.025097, 0.007244050331105379), (105.02596, 0.0006585500172305364), (105.028019, 0.0006585500172305364), (105.029568, 0.03394583873986432), (105.030431, 0.002998248471004827), (105.03249, 0.0025804814387165973), (105.034039, 0.17260036575981677), (106.025988, 0.0014725626051111385), (106.028451, 0.0016131115123158852), (106.030459, 0.005856044505252863), (106.032922, 0.006267967771419289), (106.033785, 0.0008160198806132668), (106.035844, 0.0009994161637346742), (106.037393, 0.04756832954794774), (106.038256, 0.003499783170785768), (106.040315, 0.003850598646304737), (106.041864, 0.023369277618163833), (107.033813, 0.001413387893208466), (107.036276, 0.0008160198806132668), (107.038284, 0.008227331050963668), (107.040747, 0.008684748606682719), (107.04161, 0.0013906727757838806), (107.043669, 0.0011354795666168082), (107.045218, 0.006471150231264972), (107.046081, 0.0004068377315863751), (107.04814, 0.0006150808192940338), (107.049689, 0.12409168802635434), (108.041638, 0.001795350833545003), (108.044101, 0.001605810602907014), (108.046109, 0.001019996163538284), (108.048572, 0.0011710789977235128), (108.051494, 0.00010873195457514956), (108.053043, 0.03435307625692087), (108.053906, 0.00238061116188535), (108.055965, 0.0035592280137646524), (109.049463, 0.0003261958703181823), (109.051926, 0.0002174639091502991), (109.052789, 0.00010873195457514956), (109.053934, 0.00562763357458522), (109.056397, 0.00588817455137406), (109.05726, 0.0010000567086403053), (109.059319, 0.0005773830099383702), (110.057288, 0.0016330857666772905), (110.059751, 0.0010000567086403053), (114.010564, 0.12353147876104385), (115.013918, 0.03609424972916289), (115.014781, 0.002081235972271491), (115.01684, 0.0018253646740100575), (115.054775, 0.12289876332047583), (116.014809, 0.0057723101402779526), (116.017272, 0.00711658183713664), (116.018135, 0.0005772309930282555), (116.020194, 0.0005772309930282555), (116.026214, 0.049540692854491294), (116.058129, 0.03856822345510715), (116.061051, 0.003464304885058382), (117.018163, 0.0014139213966392556), (117.020626, 0.0016326557979252918), (117.029568, 0.014576851342650719), (117.030431, 0.0010098288530943718), (117.03249, 0.001036062512829171), (117.061483, 0.007916692506060535), (117.064405, 0.0010000586587595177), (118.030459, 0.002339758570119465), (118.032922, 0.00287495680381866), (118.033785, 0.0003276317338857577), (118.035844, 0.0004012652857817208), (118.064837, 0.0011547682716804466), (119.013303, 0.041727289478160214), (119.033813, 0.0006552634677715218), (119.036276, 0.00046334124152505574), (119.037139, 0.00023167062076252565), (119.049689, 0.12340100935276455), (120.016657, 0.01139302406619951), (120.01752, 0.0012283628393010713), (120.019579, 0.0007768848730605564), (120.053043, 0.0365352614530549), (120.053906, 0.0025167610595861474), (120.055965, 0.003512093238414541), (121.017548, 0.002780824691758602), (121.020011, 0.002136433442668551), (121.020874, 0.00033640101794314345), (121.022933, 0.00019422121826513446), (121.053934, 0.005686578569157463), (121.056397, 0.007047877066786547), (121.05726, 0.0010000594088043554), (121.059319, 0.0008165450880133591), (122.020902, 0.0005493405619423678), (122.023365, 0.00047574288197131263), (122.057288, 0.0015276160197955427), (122.059751, 0.0010000594088043554), (123.021793, 0.00019422121826513446), (123.060642, 0.0005773845688788126), (130.041864, 0.0387865762682763), (131.045218, 0.01219731890878538), (131.046081, 0.0007742703974732377), (131.04814, 0.001016101865069893), (131.049689, 0.12274683847523885), (132.021128, 0.12330637956227544), (132.046109, 0.0017881008534876347), (132.048572, 0.0025549619555695004), (132.049435, 0.00036499455830664993), (132.051494, 0.00025809012727482484), (132.053043, 0.03851634981334116), (132.053906, 0.0025167606820719604), (132.055965, 0.003464306963643808), (133.024482, 0.03628587662920331), (133.025345, 0.003265570662586144), (133.027404, 0.002581660288499354), (133.028953, 0.12329830952251021), (133.049463, 0.0005161802545496417), (133.051926, 0.0004828424281325612), (133.053934, 0.005568094312730065), (133.056397, 0.007874474493635356), (133.05726, 0.0005773844822711213), (133.059319, 0.0010000592587954315), (134.025373, 0.007998981842338427), (134.027836, 0.007140519651394268), (134.028699, 0.0009998727302923031), (134.030758, 0.0009998727302923031), (134.032307, 0.03643475126539189), (134.03317, 0.003162465642353035), (134.035229, 0.0027690391822014133), (134.052817, 0.00018249727915332778), (134.057288, 0.0019149677191333348), (134.059751, 0.0011547689645422406), (135.028727, 0.0022357834751192307), (135.03119, 0.0011545535799792545), (135.032053, 0.0008163926656465358), (135.033198, 0.007724879167152918), (135.035661, 0.007095020391875296), (135.036524, 0.0010000594088043554), (135.038583, 0.0008165450880133591), (135.060642, 0.0008165449655315875), (136.036552, 0.002582142290344868), (136.037415, 0.0005773845688788126), (136.039015, 0.0014142975791099387), (136.039878, 0.0005773845688788126), (142.041864, 0.12206832045883433), (143.045218, 0.04038378497849293), (143.046081, 0.00251675300593787), (143.04814, 0.003162455522446787), (143.049689, 0.12204549617351895), (144.046109, 0.005538060642523705), (144.048572, 0.009019006449853635), (144.049435, 0.0005773827212452357), (144.051494, 0.0010000562086091265), (144.053043, 0.04058976575105422), (144.053906, 0.0021603752932932625), (144.055965, 0.003464307483289969), (144.057514, 0.29905288317871853), (145.049463, 0.002000112417218253), (145.051926, 0.0016330849501297931), (145.053934, 0.005164284580689736), (145.056397, 0.008756470746006448), (145.05726, 0.0005773845688788126), (145.059319, 0.0010000594088043554), (145.060868, 0.09888964429810151), (145.061731, 0.005657190316439752), (145.06379, 0.008944803278499157), (145.065339, 0.5176676929442573), (146.036778, 0.2103041818929387), (146.052817, 0.0008165424750648966), (146.057288, 0.0015276160197955427), (146.059751, 0.0015276160197955427), (146.061759, 0.01400083198313204), (146.064222, 0.02107255915422949), (146.065085, 0.0028285951582198766), (146.067144, 0.0020001188176087094), (146.068693, 0.17227747287402762), (146.069556, 0.010955101940948515), (146.071615, 0.015875451089430455), (147.040132, 0.06612806503135865), (147.040995, 0.006026467400326974), (147.043054, 0.005516235539997182), (147.065113, 0.004899270528080153), (147.067576, 0.003464307413272874), (147.069584, 0.022451278239463277), (147.072047, 0.03690747759011122), (147.074969, 0.005477550970474257), (148.041023, 0.01362049370911649), (148.043486, 0.013799483022363349), (148.044349, 0.002215373465254832), (148.046408, 0.0014011252036210328), (148.072938, 0.006928614826545748), (148.075401, 0.0073489059406489175), (148.07586, 0.0024496352640400766), (149.044377, 0.004318558036589935), (149.04684, 0.0026212653020769435), (150.047731, 0.0017160209073123552)])
#list3 = np.array([(1, 0.3), (2, 0.4), (3, 0.5), (4, 0.6), (5, 0.7)])
#list4 = np.array([(1, 0.4), (2, 0.1), (3, 0.2)])
list5 = np.array([(51.026928, 0.06145898604186369), (51.927979, 0.04324891610353371), (52.881775, 0.05918272729957245), (58.864681, 0.05007769233040745), (60.705158, 0.04893956295926183), (63.765099, 0.07511653849561117), (79.11158, 0.03983452799009684), (88.169312, 0.052353951072698704), (88.348816, 0.047801433588116206), (90.712204, 0.05349208044384433), (90.755814, 0.051215821701553074), (90.869606, 0.05007769233040745), (91.049934, 0.10129351403196053), (91.483452, 0.06601150352644619), (92.830452, 0.05804459792842682), (97.942276, 0.042110786732388085), (104.013573, 0.06373524478415495), (105.03331, 0.5474402275210452), (106.039116, 0.09218847906279554), (106.987946, 0.047801433588116206), (112.600822, 0.04893956295926183), (112.954582, 0.04097265736124246), (114.593483, 0.05576833918613557), (115.05407, 0.13885178327976613), (117.071342, 0.20827767491964919), (118.071045, 0.051215821701553074), (126.042984, 0.04552517484582496), (127.057457, 0.06942589163988307), (127.998978, 0.042110786732388085), (128.062027, 0.05349208044384433), (129.111496, 0.06714963289759181), (130.006332, 0.04438704547467933), (132.00415, 0.04893956295926183), (137.773453, 0.04438704547467933), (144.055267, 0.11722732522799927), (145.064026, 0.7090545982237237), (146.068298, 0.09560286717623241), (147.08812, 0.06373524478415495), (148.387192, 0.04552517484582496), (149.584854, 0.05576833918613557), (149.96727, 0.04097265736124246)])
list1 = np.array([(1, 0.5), (2, 0.3), (3, 0.2)])
list2 = np.array([(7, 0.1), (2, 0.2), (3, 0.3), (4, 0.4)])
list5 = np.array([(1, 0.6), (2, 0.3), (3, 0.1), (4, 0.4), (5, 0.5)])

list1 = [(1, 0.5), (2, 0.3), (3, 0.2)]
list2 = [(7, 0.1), (2, 0.2), (3, 0.3), (4, 0.4)]
list5 = [(1, 0.6), (2, 0.3), (3, 0.1), (4, 0.4), (5, 0.5)]

lists = [list1, list2, list5]



# Ensure all lists have the same shape
# Determine the common X-axis range
x_values = np.unique(np.concatenate([lst[:, 0] for lst in lists]))

# Pad the lists based on the common X-axis range
padded_lists = []

# Pad and fill the lists with their respective X-axis values
for lst in lists:
    padded_lst = np.zeros((len(x_values), 2))
    for i, x in enumerate(x_values):
        if x in lst[:, 0]:
            index = np.where(lst[:, 0] == x)[0][0]
            padded_lst[i] = lst[index]
    for i, x in enumerate(x_values):
        if padded_lst[i, 0] == 0:
            padded_lst[i, 0] = x
    padded_lists.append(padded_lst)

lists = padded_lists

# Define the initial weights
initial_weights = np.ones(len(lists) - 1) / (len(lists) - 1)
print("initial_weights: ", initial_weights)

print("cosine_score: ", cosine_score(initial_weights, lists))

# Minimize the cosine score function
problem = minimize(cosine_score, initial_weights, args=(lists,), method='Nelder-Mead')
optimized_weights = problem.x

# Print the optimized weights and result list
for i, weight in enumerate(optimized_weights):
    print(f"Ratio for List{i+1}: {weight:.2f}")

result_list = np.sum([w * lst for w, lst in zip(optimized_weights, lists[:-1])], axis=0)
# print("Result List:")
# print(result_list)

# Compute and print the cosine score of each list compared to list 5
cosine_scores = []
for i, lst in enumerate(lists[:-1]):
    dot_product = np.dot(lst[:, 1], lists[-1][:, 1])
    magnitude = np.linalg.norm(lst[:, 1]) * np.linalg.norm(lists[-1][:, 1])
    cosine_score = dot_product / magnitude
    cosine_scores.append(cosine_score)

print("Cosine Scores:")
for i, score in enumerate(cosine_scores):
    print(f"List{i+1} vs List5: {score:.4f}")

# Compute and print the cosine score between the combined list and list 5
combined_values = np.sum([w * lst[:, 1] for w, lst in zip(optimized_weights, lists[:-1])], axis=0)
combined_dot_product = np.dot(combined_values, lists[-1][:, 1])
combined_magnitude = np.linalg.norm(combined_values) * np.linalg.norm(lists[-1][:, 1])
combined_cosine_score = combined_dot_product / combined_magnitude
print("Combined List vs List5:")
print(f"{combined_cosine_score:.4f}")

# Plot the combined values
plt.figure(figsize=(10, 6))
plt.plot(lists[-1][:, 0], lists[-1][:, 1], label='List5', color='red', marker='o')
plt.plot(lists[-1][:, 0], combined_values, label='Combined List', color='steelblue', marker='o')
plt.xlabel('X')
plt.ylabel('Values')
plt.title('Comparison of Combined List and List5')
plt.legend()
plt.show()


plt.savefig('_intensity.png')
plt.show()