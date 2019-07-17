/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Altug Emiroglu, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
#include <math.h>
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"


#include "AuxiliaryParameters.h"
#include "MathLibrary.h"

namespace EMPIRE {
using namespace std;
using namespace MathLibrary;

/********//**
 * \brief This class manages tests the MatrixVectorMath of EMPIRE
 **************************************************************************************************/
class TestMatrixVectorMath: public CppUnit::TestFixture {
private:
    double* A;
    double* B;
    int n;
    int m;
    double tol;
public:
    /***********************************************************************************************
     * \brief Set up some test vectors and matrices
     * \author Andreas Apostolatos
     ***********/
    void setUp() {
        // Assign a tolerance
        tol = 1e-13;

        // Get the sizes of the matrices
        n = 3;
        m = 27;

        // Initialize the matrices
        A = new double[n*m];
        B = new double[n*n];

        // Assign the values to the matrices
        double ATemp[] = {0000000000000000, 0000000000000000, 0000000000000000, 1.9297090468641944e-02, -7.3617664372221947e-02,
                          1.1464583416185073e-02, 3.8378454811742241e-01, -1.4641234180077731e+00, 2.2801001906919766e-01,
                          0000000000000000, 0000000000000000, 0000000000000000, 2.1485534565930742e-01, -8.1966495161681097e-01,
                          1.2764758691094841e-01, 1.8109659197576788e+00, -6.9087659347868247e+00, 1.0759119300740561e+00,
                          0000000000000000, 0000000000000000, 0000000000000000, 1.1961082846503659e+00, -4.5631075003413031e+00,
                          7.1061921103846160e-01, -3.6250111886534158e+00, 1.3829279469124934e+01, -2.1536533305088490e+00,
                          0000000000000000, 0000000000000000, 0000000000000000, -2.5538947003129771e-18, 9.7430109059787522e-18,
                          -1.5172929243126293e-18, -5.7934837759982427e-01, 2.2101919709627977e+00, -3.4419633430324637e-01,
                          0000000000000000, 0000000000000000, 0000000000000000, -2.8435267456764674e-17, 1.0847946115857058e-16,
                          -1.6893660536512676e-17, -2.9117085748274776e+00, 1.1108057194409538e+01, -1.7298735212946448e+00,
                          0000000000000000, 0000000000000000, 0000000000000000, -1.5830026884793869e-16, 6.0390948993152508e-16,
                          -9.4047682471140848e-17, 3.4910569524273001e+00, -1.3318249165372338e+01, 2.0740698555978900e+00,
                          -3.6774185675646665e-01, 1.6062869507087991e+00, -2.1193032717698823e-01, 8.5534159708889959e-03,
                          -3.7361100471959334e-02, 4.9293497922153658e-03, 0000000000000000, 0000000000000000, 0000000000000000,
                          -1.9004557708238417e+00, 8.3011418172481051e+00, -1.0952362530839634e+00, 9.5234416192407967e-02,
                          -4.1598147498760529e-01, 5.4883774069625468e-02, 0000000000000000, 0000000000000000, 0000000000000000,
                          1.6342359994519522e+00, -7.1383007184755591e+00, 9.4181329561738636e-01, 5.3017379596506020e-01,
                          -2.3157854740217849e+00, 3.0554016078172525e-01, 0000000000000000, 0000000000000000, 0000000000000000};
        double BTemp[] = {7.7246995640955596e-01, -5.8345896723941426e-01, -7.8132120999644150e-01, 2.9043311578159398e-01,
                          -2.6531040888275976e-02, -2.5408845432800031e-01, 5.6474664381651163e-01, 8.1170914582583142e-01,
                          -5.7006685940063639e-01};

        // Assign the values to the member variables
        for(int i = 0; i < n*m; i++)
            A[i] = ATemp[i];
        for(int i = 0; i < n*n; i++)
            B[i] = BTemp[i];
    }

    /***********************************************************************************************
     * \brief Delete test vectors and matrices
     * \author Andreas Apostolatos
     ***********/
    void tearDown() {
        delete[] A;
        delete[] B;
    }

    /***********************************************************************************************
     * \brief Test
     * \author Andreas Apostolatos
     ***********/
    void testMatrixProduct() {
        // Define the product matrix
        double* C = new double[n*m];

        // Compute the product matrix
        EMPIRE::MathLibrary::computeMatrixProduct(n, n, m, B, A, C);

        // Define the expected solution from Matlab
        double CExpSol[] = {2.8732451248730057e-01, -1.2550260639292934e+00, 1.6558565966486619e-01, 8.2234573171652246e-03, -2.7676413761033723e-02,
                            5.0046407076047867e-03, 6.3448803922115005e-01, -2.4205476776654740e+00, 3.7695532723148945e-01, 1.4848664024048039e+00,
                            -6.4858581690043486e+00, 8.5573131449153117e-01, 9.1560630203048401e-02, -3.0815140007245867e-01, 5.5722069131634244e-02,
                            3.0977792430653928e+00, -1.1817909699175138e+01, 1.8404198599142094e+00, -1.2768632485115432e+00, 5.5773057546777913e+00,
                            -7.3585870372251261e-01, 5.0972168273314966e-01, -1.7154911432258808e+00, 3.1020698285018866e-01, -4.8371007189204658e+00,
                            1.8453354812155840e+01, -2.8737671502690114e+00, 9.3438959974959504e-02, -4.0813896851283560e-01, 5.3849049257628338e-02,
                            3.4311898670593621e-03, -1.1887983369274793e-02, 2.0772038131545694e-03, 1.2683445759323136e-01, -4.8386861973334794e-01,
                            7.5353547286681866e-02, 4.8288386932735833e-01, -2.1092242935020988e+00, 2.7828688667009482e-01, 3.8203141873005608e-02,
                            -1.3236175578635753e-01, 2.3127753067510241e-02, 6.0321513390299442e-01, -2.3012427362600616e+00, 3.5837579928301072e-01,
                            -4.1524049910792132e-01, 1.8137597960859089e+00, -2.3930388454898174e-01, 2.1267841558121273e-01, -7.3686317732105455e-01,
                            1.2875312440800121e-01, -1.1454446690119000e+00, 4.3698277383991231e+00, -6.8051947923710032e-01, 2.0963744535131765e-01,
                            -9.1569095728678995e-01, 1.2081445602553502e-01, 6.0219480979171961e-03, -2.0277003670021224e-02, 3.6645260521095632e-03,
                            -2.5352134121918379e-01, 9.6717425040827387e-01, -1.5061941948780225e-01, 1.0833868527033630e+00, -4.7322058451979192e+00,
                            6.2435789109729556e-01, 6.7048850821471218e-02, -2.2576577746458876e-01, 4.0801125583309110e-02, -1.3407235549177294e+00,
                            5.1148092424739433e+00, -7.9653650680533084e-01, -9.3162378372703447e-01, 4.0693086720386686e+00, -5.3689654757436656e-01,
                            3.7326362859510309e-01, -1.2568470939704892e+00, 2.2714149458790492e-01, 7.8650995429477355e-01, -3.0005054873310764e+00,
                            4.6727298052134381e-01};

        // Compare the dual product to the expected one
        for(int i = 0; i < n*m; i++)
            CPPUNIT_ASSERT(fabs(C[i] - CExpSol[i]) < tol);

        // Delete dual product matrix
        delete[] C;
    }

    /***********************************************************************************************
     * \brief Test
     * \author Andreas Apostolatos
     ***********/
    void testTransposeMatrixProduct() {
        // Define the dual poroduct matrix
        double* C = new double[m*m];

        // Compute the dual product matrix
        EMPIRE::MathLibrary::computeTransposeMatrixProduct(n, m, m, A, A, C);

        // Define the expected solution from Matlab
        double CExpSol[] = {1.3523407321069364e-01, -5.9069894573733683e-01, 7.7935652019071114e-02, -3.1454490707451353e-03, 1.3739240458023229e-02,
                            -1.8127282451913817e-03, 0000000000000000, 0000000000000000, 0000000000000000, 6.9887713384630157e-01, -3.0526773050735678e+00,
                            4.0276421329609213e-01, -3.5021681037714220e-02, 1.5297379998823565e-01, -2.0183060982166486e-02, 0000000000000000,
                            0000000000000000, 0000000000000000, -6.0097698081672091e-01, 2.6250519602982219e+00, -3.4634417004826468e-01, -1.9496709613181534e-01,
                            8.5161125006642546e-01, -1.1235990603954099e-01, 0000000000000000, 0000000000000000, 0000000000000000, -5.9069894573733683e-01,
                            2.5801577680173722e+00, -3.4042091900384258e-01, 1.3739240458023229e-02, -6.0012648152228634e-02, 7.9179502467146730e-03,
                            0000000000000000, 0000000000000000, 0000000000000000, -3.0526773050735692e+00, 1.3334015777028759e+01, -1.7592637012719703e+00,
                            1.5297379998823568e-01, -6.6818561500918905e-01, 8.8159090093689357e-02, 0000000000000000, 0000000000000000, 0000000000000000,
                            2.6250519602982232e+00, -1.1466159294322535e+01, 1.5128224067542564e+00, 8.5161125006642557e-01, -3.7198159875621837e+00,
                            4.9078517318115367e-01, 0000000000000000, 0000000000000000, 0000000000000000, 7.7935652019071114e-02, -3.4042091900384258e-01,
                            4.4914463577345277e-02, -1.8127282451913813e-03, 7.9179502467146713e-03, -1.0446787142340214e-03, 0000000000000000, 0000000000000000,
                            0000000000000000, 4.0276421329609213e-01, -1.7592637012719696e+00, 2.3211377745218306e-01, -2.0183060982166486e-02,
                            8.8159090093689330e-02, -1.1631536195283629e-02, 0000000000000000, 0000000000000000, 0000000000000000, -3.4634417004826457e-01,
                            1.5128224067542555e+00, -1.9959879987983023e-01, -1.1235990603954100e-01, 4.9078517318115367e-01, -6.4753226240180622e-02,
                            0000000000000000, 0000000000000000, 0000000000000000, -3.1454490707451353e-03, 1.3739240458023229e-02, -1.8127282451913813e-03,
                            4.4553862532601069e-04, -1.7401717629477310e-03, 2.6339588260624876e-04, 7.4059251454887693e-03, -2.8253322054553269e-02,
                            4.3999299657350825e-03, -1.6255388742132806e-02, 7.1003118996264450e-02, -9.3680112590249949e-03, 4.9606626192974235e-03,
                            -1.9375211317077569e-02, 2.9326707823942046e-03, 3.4946373189191304e-02, -1.3331908127025302e-01, 2.0761969850930227e-02,
                            1.3978300297914067e-02, -6.1056855370417239e-02, 8.0557208843293517e-03, 2.7616206792944631e-02, -1.0786257471087540e-01,
                            1.6326295295142554e-02, -6.9952168857284253e-02, 2.6686485703183654e-01, -4.1559243156921298e-02, 1.3739240458023229e-02,
                            -6.0012648152228634e-02, 7.9179502467146713e-03, -1.7401717629477310e-03, 6.8154123360969565e-03, -1.0281617869484444e-03,
                            -2.8253322054553272e-02, 1.0778534638640669e-01, -1.6785565057340121e-02, 7.1003118996264464e-02, -3.1013979346618953e-01,
                            4.0919231692002241e-02, -1.9375211317077572e-02, 7.5883344987285689e-02, -1.1447635408429826e-02, -1.3331908127025305e-01,
                            5.0860721181337676e-01, -7.9206123362261405e-02, -6.1056855370417266e-02, 2.6669477034202488e-01, -3.5187181163388308e-02,
                            -1.0786257471087543e-01, 4.2244561022092658e-01, -6.3729443219867274e-02, 2.6686485703183654e-01, -1.0180792544676993e+00,
                            1.5854692805951845e-01, -1.8127282451913817e-03, 7.9179502467146730e-03, -1.0446787142340214e-03, 2.6339588260624876e-04,
                            -1.0281617869484444e-03, 1.5573516228067945e-04, 4.3999299657350843e-03, -1.6785565057340125e-02, 2.6140398833447661e-03,
                            -9.3680112590249966e-03, 4.0919231692002234e-02, -5.3988025965661709e-03, 2.9326707823942055e-03, -1.1447635408429826e-02,
                            1.7339677283214052e-03, 2.0761969850930234e-02, -7.9206123362261405e-02, 1.2334882070802701e-02, 8.0557208843293517e-03,
                            -3.5187181163388294e-02, 4.6425271730572322e-03, 1.6326295295142561e-02, -6.3729443219867274e-02, 9.6530675501569144e-03,
                            -4.1559243156921298e-02, 1.5854692805951845e-01, -2.4690738257163503e-02, 0000000000000000, 0000000000000000, 0000000000000000,
                            7.4059251454887693e-03, -2.8253322054553272e-02, 4.3999299657350843e-03, 4.8293512200124267e-01, -1.8423790769297042e+00,
                            2.8691630998910944e-01, 0000000000000000, 0000000000000000, 0000000000000000, 8.2458161744469904e-02, -3.1457474306394673e-01,
                            4.8989171460897742e-02, 2.3819143760400490e+00, -9.0869125261984252e+00, 1.4151177919125630e+00, 0000000000000000, 0000000000000000,
                            0000000000000000, 4.5904787752404602e-01, -1.7512501500297084e+00, 2.7272467279195528e-01, -3.4137614624552945e+00, 1.3023369818276336e+01,
                            -2.0281478761202596e+00, 0000000000000000, 0000000000000000, 0000000000000000, -2.8253322054553269e-02, 1.0778534638640669e-01,
                            -1.6785565057340125e-02, -1.8423790769297042e+00, 7.0286059316671814e+00, -1.0945747829714734e+00, 0000000000000000, 0000000000000000,
                            0000000000000000, -3.1457474306394678e-01, 1.2000906505823814e+00, -1.8689182124850209e-01, -9.0869125261984252e+00, 3.4666224818735259e+01,
                            -5.3986204200819330e+00, 0000000000000000, 0000000000000000, 0000000000000000, -1.7512501500297086e+00, 6.6809525501366158e+00,
                            -1.0404342281676198e+00, 1.3023369818276333e+01, -4.9683659297507845e+01, 7.7373068175268660e+00, 0000000000000000, 0000000000000000,
                            0000000000000000, 4.3999299657350825e-03, -1.6785565057340121e-02, 2.6140398833447661e-03, 2.8691630998910944e-01, -1.0945747829714734e+00,
                            1.7045968534372802e-01, 0000000000000000, 0000000000000000, 0000000000000000, 4.8989171460897742e-02, -1.8689182124850207e-01,
                            2.9104928725702417e-02, 1.4151177919125630e+00, -5.3986204200819330e+00, 8.4073482453082837e-01, 0000000000000000, 0000000000000000,
                            0000000000000000, 2.7272467279195528e-01, -1.0404342281676195e+00, 1.6202829985981784e-01, -2.0281478761202596e+00, 7.7373068175268696e+00,
                            -1.2049417783434211e+00, 6.9887713384630157e-01, -3.0526773050735692e+00, 4.0276421329609213e-01, -1.6255388742132806e-02, 7.1003118996264464e-02,
                            -9.3680112590249966e-03, 0000000000000000, 0000000000000000, 0000000000000000, 3.6117321368576421e+00, -1.5775952871016273e+01, 2.0814480575888998e+00,
                            -1.8098879583390123e-01, 7.9055439469600797e-01, -1.0430418515521164e-01, 0000000000000000, 0000000000000000, 0000000000000000, -3.1057932360465310e+00,
                            1.3566024794302852e+01, -1.7898745126946827e+00, -1.0075718500813806e+00, 4.4010478680947269e+00, -5.8066556177607420e-01, 0000000000000000,
                            0000000000000000, 0000000000000000, -3.0526773050735678e+00, 1.3334015777028759e+01, -1.7592637012719696e+00, 7.1003118996264450e-02,
                            -3.1013979346618953e-01, 4.0919231692002234e-02, 0000000000000000, 0000000000000000, 0000000000000000, -1.5775952871016273e+01, 6.8908955470065166e+01,
                            -9.0917114602414184e+00, 7.9055439469600786e-01, -3.4531212172201569e+00, 4.5559799201776519e-01, 0000000000000000, 0000000000000000, 0000000000000000,
                            1.3566024794302852e+01, -5.9256046598229659e+01, 7.8181257322897375e+00, 4.4010478680947260e+00, -1.9223663638177964e+01, 2.5363322055138888e+00,
                            0000000000000000, 0000000000000000, 0000000000000000, 4.0276421329609213e-01, -1.7592637012719703e+00, 2.3211377745218306e-01, -9.3680112590249949e-03,
                            4.0919231692002241e-02, -5.3988025965661709e-03, 0000000000000000, 0000000000000000, 0000000000000000, 2.0814480575888998e+00, -9.0917114602414184e+00,
                            1.1995424500693996e+00, -1.0430418515521164e-01, 4.5559799201776530e-01, -6.0110699067123387e-02, 0000000000000000, 0000000000000000, 0000000000000000,
                            -1.7898745126946822e+00, 7.8181257322897357e+00, -1.0315080649966455e+00, -5.8066556177607431e-01, 2.5363322055138897e+00, -3.3463866086124849e-01,
                            0000000000000000, 0000000000000000, 0000000000000000, -3.5021681037714220e-02, 1.5297379998823568e-01, -2.0183060982166486e-02, 4.9606626192974235e-03,
                            -1.9375211317077572e-02, 2.9326707823942055e-03, 8.2458161744469904e-02, -3.1457474306394678e-01, 4.8989171460897742e-02, -1.8098879583390123e-01,
                            7.9055439469600786e-01, -1.0430418515521164e-01, 5.5232413585889248e-02, -2.1572514942175075e-01, 3.2652590590285108e-02, 3.8909570866676174e-01,
                            -1.4843852929978716e+00, 2.3116542963503395e-01, 1.5563551132842302e-01, -6.7981190152986626e-01, 8.9693039370369529e-02, 3.0748105088376099e-01,
                            -1.2009505169117323e+00, 1.8177827505519506e-01, -7.7885303195698663e-01, 2.9712946205580018e+00, -4.6272393075679746e-01, 1.5297379998823565e-01,
                            -6.6818561500918905e-01, 8.8159090093689330e-02, -1.9375211317077569e-02, 7.5883344987285689e-02, -1.1447635408429826e-02, -3.1457474306394673e-01,
                            1.2000906505823814e+00, -1.8689182124850207e-01, 7.9055439469600797e-01, -3.4531212172201569e+00, 4.5559799201776530e-01, -2.1572514942175075e-01,
                            8.4489122044185272e-01, -1.2745888643973449e-01, -1.4843852929978716e+00, 5.6628732956689154e+00, -8.8188730010810112e-01, -6.7981190152986637e-01,
                            2.9694008617765455e+00, -3.9177688387385790e-01, -1.2009505169117325e+00, 4.7035431457280144e+00, -7.0956870798374916e-01, 2.9712946205580009e+00,
                            -1.1335375686955649e+01, 1.7652741529509197e+00, -2.0183060982166486e-02, 8.8159090093689357e-02, -1.1631536195283629e-02, 2.9326707823942046e-03,
                            -1.1447635408429826e-02, 1.7339677283214052e-03, 4.8989171460897742e-02, -1.8689182124850209e-01, 2.9104928725702417e-02, -1.0430418515521164e-01,
                            4.5559799201776519e-01, -6.0110699067123387e-02, 3.2652590590285108e-02, -1.2745888643973449e-01, 1.9306135100313822e-02, 2.3116542963503398e-01,
                            -8.8188730010810112e-01, 1.3733756160265434e-01, 8.9693039370369515e-02, -3.9177688387385773e-01, 5.1690268132434013e-02, 1.8177827505519509e-01,
                            -7.0956870798374916e-01, 1.0747802465516286e-01, -4.6272393075679735e-01, 1.7652741529509197e+00, -2.7490865068218184e-01, 0000000000000000,
                            0000000000000000, 0000000000000000, 3.4946373189191304e-02, -1.3331908127025305e-01, 2.0761969850930234e-02, 2.3819143760400490e+00, -9.0869125261984252e+00,
                            1.4151177919125630e+00, 0000000000000000, 0000000000000000, 0000000000000000, 3.8909570866676174e-01, -1.4843852929978716e+00, 2.3116542963503398e-01,
                            1.1757644387247636e+01, -4.4854965038118046e+01, 6.9853274033454431e+00, 0000000000000000, 0000000000000000, 0000000000000000, 2.1661113398416298e+00,
                            -8.2636321713087533e+00, 1.2869071731157440e+00, -1.6729712184985260e+01, 6.3823214109893357e+01, -9.9392797678597109e+00, 0000000000000000, 0000000000000000,
                            0000000000000000, -1.3331908127025302e-01, 5.0860721181337676e-01, -7.9206123362261405e-02, -9.0869125261984252e+00, 3.4666224818735259e+01,
                            -5.3986204200819330e+00, 0000000000000000, 0000000000000000, 0000000000000000, -1.4843852929978716e+00, 5.6628732956689154e+00, -8.8188730010810112e-01,
                            -4.4854965038118046e+01, 1.7111998137594435e+02, -2.6648757704961923e+01, 0000000000000000, 0000000000000000, 0000000000000000, -8.2636321713087515e+00,
                            3.1525441655128262e+01, -4.9095017978276134e+00, 6.3823214109893328e+01, -2.4348312835729018e+02, 3.7917973346341824e+01, 0000000000000000, 0000000000000000,
                            0000000000000000, 2.0761969850930227e-02, -7.9206123362261405e-02, 1.2334882070802701e-02, 1.4151177919125630e+00, -5.3986204200819330e+00, 8.4073482453082837e-01,
                            0000000000000000, 0000000000000000, 0000000000000000, 2.3116542963503395e-01, -8.8188730010810112e-01, 1.3733756160265434e-01, 6.9853274033454431e+00, -2.6648757704961923e+01,
                            4.1500488809520144e+00, 0000000000000000, 0000000000000000, 0000000000000000, 1.2869071731157438e+00, -4.9095017978276134e+00, 7.6456368689609433e-01, -9.9392797678597091e+00,
                            3.7917973346341832e+01, -5.9050198360523920e+00, -6.0097698081672091e-01, 2.6250519602982232e+00, -3.4634417004826457e-01, 1.3978300297914067e-02, -6.1056855370417266e-02,
                            8.0557208843293517e-03, 0000000000000000, 0000000000000000, 0000000000000000, -3.1057932360465310e+00, 1.3566024794302852e+01, -1.7898745126946822e+00,
                            1.5563551132842302e-01, -6.7981190152986637e-01, 8.9693039370369515e-02, 0000000000000000, 0000000000000000, 0000000000000000, 2.6707273019047211e+00,
                            -1.1665668009046493e+01, 1.5391451924604163e+00, 8.6642910333219558e-01, -3.7845399886543047e+00, 4.9932473002783290e-01, 0000000000000000, 0000000000000000,
                            0000000000000000, 2.6250519602982219e+00, -1.1466159294322535e+01, 1.5128224067542555e+00, -6.1056855370417239e-02, 2.6669477034202488e-01, -3.5187181163388294e-02,
                            0000000000000000, 0000000000000000, 0000000000000000, 1.3566024794302852e+01, -5.9256046598229659e+01, 7.8181257322897357e+00, -6.7981190152986626e-01,
                            2.9694008617765455e+00, -3.9177688387385773e-01, 0000000000000000, 0000000000000000, 0000000000000000, -1.1665668009046493e+01, 5.0955337147388683e+01,
                            -6.7229465247754234e+00, -3.7845399886543039e+00, 1.6530773113044969e+01, -2.1810375492313270e+00, 0000000000000000, 0000000000000000, 0000000000000000,
                            -3.4634417004826468e-01, 1.5128224067542564e+00, -1.9959879987983023e-01, 8.0557208843293517e-03, -3.5187181163388308e-02, 4.6425271730572322e-03,
                            0000000000000000, 0000000000000000, 0000000000000000, -1.7898745126946827e+00, 7.8181257322897375e+00, -1.0315080649966455e+00, 8.9693039370369529e-02,
                            -3.9177688387385790e-01, 5.1690268132434013e-02, 0000000000000000, 0000000000000000, 0000000000000000, 1.5391451924604163e+00, -6.7229465247754234e+00,
                            8.8701228380168240e-01, 4.9932473002783312e-01, -2.1810375492313283e+00, 2.8776178576930278e-01, 0000000000000000, 0000000000000000, 0000000000000000,
                            -1.9496709613181534e-01, 8.5161125006642557e-01, -1.1235990603954100e-01, 2.7616206792944631e-02, -1.0786257471087543e-01, 1.6326295295142561e-02,
                            4.5904787752404602e-01, -1.7512501500297086e+00, 2.7272467279195528e-01, -1.0075718500813806e+00, 4.4010478680947260e+00, -5.8066556177607431e-01,
                            3.0748105088376099e-01, -1.2009505169117325e+00, 1.8177827505519509e-01, 2.1661113398416298e+00, -8.2636321713087515e+00, 1.2869071731157438e+00,
                            8.6642910333219558e-01, -3.7845399886543039e+00, 4.9932473002783312e-01, 1.7117592825372419e+00, -6.6857394603113312e+00, 1.0119669124162327e+00,
                            -4.3359059146986221e+00, 1.6541315743765551e+01, -2.5760025908864872e+00, 8.5161125006642546e-01, -3.7198159875621837e+00, 4.9078517318115367e-01,
                            -1.0786257471087540e-01, 4.2244561022092658e-01, -6.3729443219867274e-02, -1.7512501500297084e+00, 6.6809525501366158e+00, -1.0404342281676195e+00,
                            4.4010478680947269e+00, -1.9223663638177964e+01, 2.5363322055138897e+00, -1.2009505169117323e+00, 4.7035431457280144e+00, -7.0956870798374916e-01,
                            -8.2636321713087533e+00, 3.1525441655128262e+01, -4.9095017978276134e+00, -3.7845399886543047e+00, 1.6530773113044969e+01, -2.1810375492313283e+00,
                            -6.6857394603113312e+00, 2.6184812421361361e+01, -3.9501973178448235e+00, 1.6541315743765548e+01, -6.3104488869879987e+01, 9.8273516655799575e+00,
                            -1.1235990603954099e-01, 4.9078517318115367e-01, -6.4753226240180622e-02, 1.6326295295142554e-02, -6.3729443219867274e-02, 9.6530675501569144e-03,
                            2.7272467279195528e-01, -1.0404342281676198e+00, 1.6202829985981784e-01, -5.8066556177607420e-01, 2.5363322055138888e+00, -3.3463866086124849e-01,
                            1.8177827505519506e-01, -7.0956870798374916e-01, 1.0747802465516286e-01, 1.2869071731157440e+00, -4.9095017978276134e+00, 7.6456368689609433e-01,
                            4.9932473002783290e-01, -2.1810375492313270e+00, 2.8776178576930278e-01, 1.0119669124162327e+00, -3.9501973178448235e+00, 5.9833445294744814e-01,
                            -2.5760025908864868e+00, 9.8273516655799575e+00, -1.5304274305765537e+00, 0000000000000000, 0000000000000000, 0000000000000000, -6.9952168857284253e-02,
                            2.6686485703183654e-01, -4.1559243156921298e-02, -3.4137614624552945e+00, 1.3023369818276333e+01, -2.0281478761202596e+00, 0000000000000000,
                            0000000000000000, 0000000000000000, -7.7885303195698663e-01, 2.9712946205580009e+00, -4.6272393075679735e-01, -1.6729712184985260e+01, 6.3823214109893328e+01,
                            -9.9392797678597091e+00, 0000000000000000, 0000000000000000, 0000000000000000, -4.3359059146986221e+00, 1.6541315743765548e+01, -2.5760025908864868e+00,
                            2.5328184762953438e+01, -9.6626059149525048e+01, 1.5047713408780172e+01, 0000000000000000, 0000000000000000, 0000000000000000, 2.6686485703183654e-01,
                            -1.0180792544676993e+00, 1.5854692805951845e-01, 1.3023369818276336e+01, -4.9683659297507845e+01, 7.7373068175268696e+00, 0000000000000000, 0000000000000000,
                            0000000000000000, 2.9712946205580018e+00, -1.1335375686955649e+01, 1.7652741529509197e+00, 6.3823214109893357e+01, -2.4348312835729018e+02, 3.7917973346341832e+01,
                            0000000000000000, 0000000000000000, 0000000000000000, 1.6541315743765551e+01, -6.3104488869879987e+01, 9.8273516655799575e+00, -9.6626059149525048e+01, 3.6862473146610137e+02,
                            -5.7406452910459080e+01, 0000000000000000, 0000000000000000, 0000000000000000, -4.1559243156921298e-02, 1.5854692805951845e-01, -2.4690738257163503e-02,
                            -2.0281478761202596e+00, 7.7373068175268660e+00, -1.2049417783434211e+00, 0000000000000000, 0000000000000000, 0000000000000000, -4.6272393075679746e-01,
                            1.7652741529509197e+00, -2.7490865068218184e-01, -9.9392797678597109e+00, 3.7917973346341824e+01, -5.9050198360523920e+00, 0000000000000000, 0000000000000000, 0000000000000000,
                            -2.5760025908864872e+00, 9.8273516655799575e+00, -1.5304274305765537e+00, 1.5047713408780172e+01, -5.7406452910459080e+01, 8.9399884339117097e+00};

        // Compare the dual product to the expected one
        for(int i = 0; i < m*m; i++)
            CPPUNIT_ASSERT(fabs(C[i] - CExpSol[i]) < tol);

        // Delete dual product matrix
        delete[] C;
    }


    void testMatrixProducts4Leakage() {

        double* C = new double[n*m];
        double* D = new double[m*m];
        for(int i = 0; i < 1e10; i++){
            EMPIRE::MathLibrary::computeMatrixProduct(n, n, m, B, A, C);
            EMPIRE::MathLibrary::computeTransposeMatrixProduct(n, m, m, A, A, D);
        }
        delete[] C;
        delete[] D;

    }


    CPPUNIT_TEST_SUITE(TestMatrixVectorMath);
    CPPUNIT_TEST(testMatrixProduct);
    CPPUNIT_TEST(testTransposeMatrixProduct);
//    CPPUNIT_TEST(testMatrixProducts4Leakage);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMatrixVectorMath);
