// Class: ReadBDT
// Automatically generated by MethodBase::MakeClass
//

/* configuration options =====================================================

#GEN -*-*-*-*-*-*-*-*-*-*-*- general info -*-*-*-*-*-*-*-*-*-*-*-

Method         : BDT::BDT
TMVA Release   : 4.1.4         [262404]
ROOT Release   : 5.34/04       [336388]
Creator        : mukul
Date           : Tue Apr 23 14:51:52 2013
Host           : Darwin proof.cern.ch 12.2.0 Darwin Kernel Version 12.2.0: Sat Aug 25 00:48:52 PDT 2012; root:xnu-2050.18.24~1/RELEASE_X86_64 x86_64
Dir            : /Users/mukul/Dropbox/Rice/2012 Fall/Thesis/squark-search
Training events: 276
Analysis type  : [Classification]


#OPT -*-*-*-*-*-*-*-*-*-*-*-*- options -*-*-*-*-*-*-*-*-*-*-*-*-

# Set by User:
V: "False" [Verbose output (short form of "VerbosityLevel" below - overrides the latter one)]
H: "True" [Print method-specific help message]
NTrees: "70" [Number of trees in the forest]
AdaBoostBeta: "5.000000e-01" [Parameter for AdaBoost algorithm]
nEventsMin: "100" [Minimum number of events required in a leaf node (default: Classification: max(40, N_train/(Nvar^2)/10), Regression: 10)]
PruneStrength: "5.000000e-01" [Pruning strength]
# Default:
VerbosityLevel: "Default" [Verbosity level]
VarTransform: "None" [List of variable transformations performed before training, e.g., "D_Background,P_Signal,G,N_AllClasses" for: "Decorrelation, PCA-transformation, Gaussianisation, Normalisation, each for the given class of events ('AllClasses' denotes all events of all classes, if no class indication is given, 'All' is assumed)"]
CreateMVAPdfs: "False" [Create PDFs for classifier outputs (signal and background)]
IgnoreNegWeightsInTraining: "False" [Events with negative weights are ignored in the training (but are included for testing and performance evaluation)]
RenormByClass: "False" [Individually re-normalize each event class to the original size after boosting]
BoostType: "AdaBoost" [Boosting type for the trees in the forest]
AdaBoostR2Loss: "quadratic" [Type of Loss function in AdaBoostR2t (Linear,Quadratic or Exponential)]
UseBaggedGrad: "False" [Use only a random subsample of all events for growing the trees in each iteration. (Only valid for GradBoost)]
GradBaggingFraction: "6.000000e-01" [Defines the fraction of events to be used in each iteration when UseBaggedGrad=kTRUE. (Only valid for GradBoost)]
Shrinkage: "1.000000e+00" [Learning rate for GradBoost algorithm]
UseRandomisedTrees: "False" [Choose at each node splitting a random set of variables]
UseNvars: "4" [Number of variables used if randomised tree option is chosen]
UsePoissonNvars: "True" [Interpret "UseNvars" not as fixed number but as mean of a Possion distribution in each split]
UseNTrainEvents: "276" [Number of randomly picked training events used in randomised (and bagged) trees]
UseWeightedTrees: "True" [Use weighted trees or simple average in classification from the forest]
UseYesNoLeaf: "True" [Use Sig or Bkg categories, or the purity=S/(S+B) as classification of the leaf node]
NodePurityLimit: "5.000000e-01" [In boosting/pruning, nodes with purity > NodePurityLimit are signal; background otherwise.]
SeparationType: "giniindex" [Separation criterion for node splitting]
nCuts: "20" [Number of steps during node cut optimisation]
UseFisherCuts: "False" [Use multivariate splits using the Fisher criterion]
MinLinCorrForFisher: "8.000000e-01" [The minimum linear correlation between two variables demanded for use in Fisher criterion in node splitting]
UseExclusiveVars: "False" [Variables already used in fisher criterion are not anymore analysed individually for node splitting]
PruneMethod: "nopruning" [Method used for pruning (removal) of statistically insignificant branches]
PruneBeforeBoost: "False" [Flag to prune the tree before applying boosting algorithm]
PruningValFraction: "5.000000e-01" [Fraction of events to use for optimizing automatic pruning.]
NNodesMax: "100000" [Max number of nodes in tree]
MaxDepth: "3" [Max depth of the decision tree allowed]
DoBoostMonitor: "False" [Create control plot with ROC integral vs tree number]
NegWeightTreatment: "inverseboostnegweights" [How to treat events with negative weights in the BDT training (particular the boosting) : Ignore;  Boost With inverse boostweight; Pair events with negative and positive weights in traning sample and *annihilate* them (experimental!); Randomly pair events with negative and positive weights in leaf node and do not boost them (experimental!) ]
##


#VAR -*-*-*-*-*-*-*-*-*-*-*-* variables *-*-*-*-*-*-*-*-*-*-*-*-

NVar 12
m3a                           m3a                           m3a                           m3a                                                             'F'    [125.568283081,1638.4765625]
m2a                           m2a                           m2a                           m2a                                                             'F'    [27.9763832092,1532.1940918]
angles_b_a                    angles_b_a                    angles_b_a                    angles_b_a                                                      'F'    [0.00122600002214,3.13862609863]
angles_j1_a                   angles_j1_a                   angles_j1_a                   angles_j1_a                                                     'F'    [0.00630699982867,3.13826608658]
angles_j2_a                   angles_j2_a                   angles_j2_a                   angles_j2_a                                                     'F'    [0.0762730017304,3.14012908936]
m3b                           m3b                           m3b                           m3b                                                             'F'    [95.4901199341,1796.98071289]
m2b                           m2b                           m2b                           m2b                                                             'F'    [18.6513938904,1735.52258301]
angles_b_b                    angles_b_b                    angles_b_b                    angles_b_b                                                      'F'    [0.00222199992277,3.13741993904]
angles_j1_b                   angles_j1_b                   angles_j1_b                   angles_j1_b                                                     'F'    [0.00158100004774,3.08318305016]
angles_j2_b                   angles_j2_b                   angles_j2_b                   angles_j2_b                                                     'F'    [0.000953999988269,3.13477706909]
mT_b                          mT_b                          mT_b                          mT_b                                                            'F'    [0.0234629996121,1683126.875]
eT_missing                    eT_missing                    eT_missing                    eT_missing                                                      'F'    [100.249580383,1035.71020508]
NSpec 0


============================================================================ */

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#define NN new BDTNode
   
#ifndef BDTNode__def
#define BDTNode__def
   
class BDTNode {
   
public:
   
   // constructor of an essentially "empty" node floating in space
   BDTNode ( BDTNode* left,BDTNode* right,
                          int selector, double cutValue, bool cutType, 
                          int nodeType, double purity, double response ) :
   fLeft         ( left         ),
   fRight        ( right        ),
   fSelector     ( selector     ),
   fCutValue     ( cutValue     ),
   fCutType      ( cutType      ),
   fNodeType     ( nodeType     ),
   fPurity       ( purity       ),
   fResponse     ( response     ){
   }

   virtual ~BDTNode();

   // test event if it decends the tree at this node to the right
   virtual bool GoesRight( const std::vector<double>& inputValues ) const;
   BDTNode* GetRight( void )  {return fRight; };

   // test event if it decends the tree at this node to the left 
   virtual bool GoesLeft ( const std::vector<double>& inputValues ) const;
   BDTNode* GetLeft( void ) { return fLeft; };   

   // return  S/(S+B) (purity) at this node (from  training)

   double GetPurity( void ) const { return fPurity; } 
   // return the node type
   int    GetNodeType( void ) const { return fNodeType; }
   double GetResponse(void) const {return fResponse;}

private:

   BDTNode*   fLeft;     // pointer to the left daughter node
   BDTNode*   fRight;    // pointer to the right daughter node
   int                     fSelector; // index of variable used in node selection (decision tree)   
   double                  fCutValue; // cut value appplied on this node to discriminate bkg against sig
   bool                    fCutType;  // true: if event variable > cutValue ==> signal , false otherwise
   int                     fNodeType; // Type of node: -1 == Bkg-leaf, 1 == Signal-leaf, 0 = internal 
   double                  fPurity;   // Purity of node from training
   double                  fResponse; // Regression response value of node
}; 
   
//_______________________________________________________________________
   BDTNode::~BDTNode()
{
   if (fLeft  != NULL) delete fLeft;
   if (fRight != NULL) delete fRight;
}; 
   
//_______________________________________________________________________
bool BDTNode::GoesRight( const std::vector<double>& inputValues ) const
{
   // test event if it decends the tree at this node to the right
   bool result;
     result = (inputValues[fSelector] > fCutValue );
   if (fCutType == true) return result; //the cuts are selecting Signal ;
   else return !result;
}
   
//_______________________________________________________________________
bool BDTNode::GoesLeft( const std::vector<double>& inputValues ) const
{
   // test event if it decends the tree at this node to the left
   if (!this->GoesRight(inputValues)) return true;
   else return false;
}
   
#endif
   
#ifndef IClassifierReader__def
#define IClassifierReader__def

class IClassifierReader {

 public:

   // constructor
   IClassifierReader() : fStatusIsClean( true ) {}
   virtual ~IClassifierReader() {}

   // return classifier response
   virtual double GetMvaValue( const std::vector<double>& inputValues ) const = 0;

   // returns classifier status
   bool IsStatusClean() const { return fStatusIsClean; }

 protected:

   bool fStatusIsClean;
};

#endif

class ReadBDT : public IClassifierReader {

 public:

   // constructor
   ReadBDT( std::vector<std::string>& theInputVars ) 
      : IClassifierReader(),
        fClassName( "ReadBDT" ),
        fNvars( 12 ),
        fIsNormalised( false )
   {      
      // the training input variables
      const char* inputVars[] = { "m3a", "m2a", "angles_b_a", "angles_j1_a", "angles_j2_a", "m3b", "m2b", "angles_b_b", "angles_j1_b", "angles_j2_b", "mT_b", "eT_missing" };

      // sanity checks
      if (theInputVars.size() <= 0) {
         std::cout << "Problem in class \"" << fClassName << "\": empty input vector" << std::endl;
         fStatusIsClean = false;
      }

      if (theInputVars.size() != fNvars) {
         std::cout << "Problem in class \"" << fClassName << "\": mismatch in number of input values: "
                   << theInputVars.size() << " != " << fNvars << std::endl;
         fStatusIsClean = false;
      }

      // validate input variables
      for (size_t ivar = 0; ivar < theInputVars.size(); ivar++) {
         if (theInputVars[ivar] != inputVars[ivar]) {
            std::cout << "Problem in class \"" << fClassName << "\": mismatch in input variable names" << std::endl
                      << " for variable [" << ivar << "]: " << theInputVars[ivar].c_str() << " != " << inputVars[ivar] << std::endl;
            fStatusIsClean = false;
         }
      }

      // initialize min and max vectors (for normalisation)
      fVmin[0] = 125.568283081055;
      fVmax[0] = 1638.4765625;
      fVmin[1] = 27.9763832092285;
      fVmax[1] = 1532.19409179688;
      fVmin[2] = 0.00122600002214313;
      fVmax[2] = 3.13862609863281;
      fVmin[3] = 0.00630699982866645;
      fVmax[3] = 3.13826608657837;
      fVmin[4] = 0.076273001730442;
      fVmax[4] = 3.14012908935547;
      fVmin[5] = 95.490119934082;
      fVmax[5] = 1796.98071289062;
      fVmin[6] = 18.6513938903809;
      fVmax[6] = 1735.52258300781;
      fVmin[7] = 0.00222199992276728;
      fVmax[7] = 3.13741993904114;
      fVmin[8] = 0.00158100004773587;
      fVmax[8] = 3.08318305015564;
      fVmin[9] = 0.000953999988269061;
      fVmax[9] = 3.1347770690918;
      fVmin[10] = 0.023462999612093;
      fVmax[10] = 1683126.875;
      fVmin[11] = 100.249580383301;
      fVmax[11] = 1035.71020507812;

      // initialize input variable types
      fType[0] = 'F';
      fType[1] = 'F';
      fType[2] = 'F';
      fType[3] = 'F';
      fType[4] = 'F';
      fType[5] = 'F';
      fType[6] = 'F';
      fType[7] = 'F';
      fType[8] = 'F';
      fType[9] = 'F';
      fType[10] = 'F';
      fType[11] = 'F';

      // initialize constants
      Initialize();

   }

   // destructor
   virtual ~ReadBDT() {
      Clear(); // method-specific
   }

   // the classifier response
   // "inputValues" is a vector of input values in the same order as the 
   // variables given to the constructor
   double GetMvaValue( const std::vector<double>& inputValues ) const;

 private:

   // method-specific destructor
   void Clear();

   // common member variables
   const char* fClassName;

   const size_t fNvars;
   size_t GetNvar()           const { return fNvars; }
   char   GetType( int ivar ) const { return fType[ivar]; }

   // normalisation of input variables
   const bool fIsNormalised;
   bool IsNormalised() const { return fIsNormalised; }
   double fVmin[12];
   double fVmax[12];
   double NormVariable( double x, double xmin, double xmax ) const {
      // normalise to output range: [-1, 1]
      return 2*(x - xmin)/(xmax - xmin) - 1.0;
   }

   // type of input variable: 'F' or 'I'
   char   fType[12];

   // initialize internal variables
   void Initialize();
   double GetMvaValue__( const std::vector<double>& inputValues ) const;

   // private members (method specific)
   std::vector<BDTNode*> fForest;       // i.e. root nodes of decision trees
   std::vector<double>                fBoostWeights; // the weights applied in the individual boosts
};

double ReadBDT::GetMvaValue__( const std::vector<double>& inputValues ) const
{
   double myMVA = 0;
   double norm  = 0;
   for (unsigned int itree=0; itree<fForest.size(); itree++){
      BDTNode *current = fForest[itree];
      while (current->GetNodeType() == 0) { //intermediate node
         if (current->GoesRight(inputValues)) current=(BDTNode*)current->GetRight();
         else current=(BDTNode*)current->GetLeft();
      }
      myMVA += fBoostWeights[itree] *  current->GetNodeType();
      norm  += fBoostWeights[itree];
   }
   return myMVA /= norm;
};

void ReadBDT::Initialize()
{
  // itree = 0
  fBoostWeights.push_back(0.379552574175871);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.654321,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.280702,-99) , 
11, 233.887, 0, 0, 0.5,-99)    );
  // itree = 1
  fBoostWeights.push_back(0.209524);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.571368,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.34875,-99) , 
7, 1.19658, 0, 0, 0.482504,-99)    );
  // itree = 2
  fBoostWeights.push_back(0.226328);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.595784,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.376982,-99) , 
11, 278.433, 0, 0, 0.471363,-99)    );
  // itree = 3
  fBoostWeights.push_back(0.157963);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.550273,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.385103,-99) , 
5, 257.537, 1, 0, 0.478561,-99)    );
  // itree = 4
  fBoostWeights.push_back(0.138932);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.535205,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.380801,-99) , 
8, 1.02878, 0, 0, 0.47304,-99)    );
  // itree = 5
  fBoostWeights.push_back(0.121145);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.522213,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.384813,-99) , 
11, 233.887, 0, 0, 0.46596,-99)    );
  // itree = 6
  fBoostWeights.push_back(0.126736);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.521469,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.388028,-99) , 
7, 1.34588, 0, 0, 0.460189,-99)    );
  // itree = 7
  fBoostWeights.push_back(0.114516);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.511902,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.375394,-99) , 
9, 1.19479, 0, 0, 0.457268,-99)    );
  // itree = 8
  fBoostWeights.push_back(0.098465);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.500297,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.394818,-99) , 
6, 100.407, 1, 0, 0.451243,-99)    );
  // itree = 9
  fBoostWeights.push_back(0.111318);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.505523,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.407472,-99) , 
11, 278.433, 0, 0, 0.449279,-99)    );
  // itree = 10
  fBoostWeights.push_back(0.0941101);
  fForest.push_back( 
NN(
0, 
0, 
-1, 1.19479, 0, -1, 0.453083,-99)    );
  // itree = 11
  fBoostWeights.push_back(0.0760047);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.511855,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.421602,-99) , 
9, 1.19479, 0, 0, 0.47649,-99)    );
  // itree = 12
  fBoostWeights.push_back(0.0784854);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.510621,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.427377,-99) , 
7, 1.34588, 0, 0, 0.472302,-99)    );
  // itree = 13
  fBoostWeights.push_back(0.0661888);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503013,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.419117,-99) , 
1, 171.235, 1, 0, 0.470656,-99)    );
  // itree = 14
  fBoostWeights.push_back(0.0811856);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.508257,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.433919,-99) , 
8, 1.61575, 0, 0, 0.466801,-99)    );
  // itree = 15
  fBoostWeights.push_back(0.062001);
  fForest.push_back( 
NN(
0, 
0, 
-1, 257.537, 1, -1, 0.469039,-99)    );
  // itree = 16
  fBoostWeights.push_back(0.058784);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.51233,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.44878,-99) , 
5, 257.537, 1, 0, 0.484505,-99)    );
  // itree = 17
  fBoostWeights.push_back(0.0558464);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.50956,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.449572,-99) , 
2, 1.64463, 0, 0, 0.482648,-99)    );
  // itree = 18
  fBoostWeights.push_back(0.0507364);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.505471,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.445231,-99) , 
11, 233.887, 0, 0, 0.481185,-99)    );
  // itree = 19
  fBoostWeights.push_back(0.0729625);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.520974,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.454882,-99) , 
4, 2.11884, 1, 0, 0.478701,-99)    );
  // itree = 20
  fBoostWeights.push_back(0.0523185);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.509655,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.456583,-99) , 
0, 269.655, 1, 0, 0.483749,-99)    );
  // itree = 21
  fBoostWeights.push_back(0.0457637);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.505881,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.457664,-99) , 
6, 100.407, 1, 0, 0.483416,-99)    );
  // itree = 22
  fBoostWeights.push_back(0.0381372);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.50137,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.453019,-99) , 
9, 1.19479, 0, 0, 0.482618,-99)    );
  // itree = 23
  fBoostWeights.push_back(0.0535563);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.508432,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.45954,-99) , 
11, 278.433, 0, 0, 0.480465,-99)    );
  // itree = 24
  fBoostWeights.push_back(0.0403831);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502197,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.454994,-99) , 
7, 1.19658, 0, 0, 0.482368,-99)    );
  // itree = 25
  fBoostWeights.push_back(0.0385416);
  fForest.push_back( 
NN(
0, 
0, 
-1, 1.02878, 0, -1, 0.480739,-99)    );
  // itree = 26
  fBoostWeights.push_back(0.0332046);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.505898,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.467998,-99) , 
8, 1.02878, 0, 0, 0.490366,-99)    );
  // itree = 27
  fBoostWeights.push_back(0.0423873);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.511266,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.470845,-99) , 
3, 2.69084, 0, 0, 0.488863,-99)    );
  // itree = 28
  fBoostWeights.push_back(0.0490725);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.517217,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.470161,-99) , 
3, 2.39256, 1, 0, 0.490003,-99)    );
  // itree = 29
  fBoostWeights.push_back(0.0340626);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.50731,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.467688,-99) , 
8, 1.7625, 1, 0, 0.491915,-99)    );
  // itree = 30
  fBoostWeights.push_back(0.0283945);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503428,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.468731,-99) , 
9, 1.19479, 0, 0, 0.490011,-99)    );
  // itree = 31
  fBoostWeights.push_back(0.0377002);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.508084,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.472438,-99) , 
3, 2.69084, 0, 0, 0.488398,-99)    );
  // itree = 32
  fBoostWeights.push_back(0.0424468);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.512561,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.47249,-99) , 
3, 2.39256, 1, 0, 0.489375,-99)    );
  // itree = 33
  fBoostWeights.push_back(0.0345659);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.50927,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476219,-99) , 
3, 2.69084, 0, 0, 0.491036,-99)    );
  // itree = 34
  fBoostWeights.push_back(0.0307269);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.506168,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.471398,-99) , 
8, 1.02878, 0, 0, 0.491924,-99)    );
  // itree = 35
  fBoostWeights.push_back(0.0259553);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502865,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.471075,-99) , 
8, 1.7625, 1, 0, 0.490532,-99)    );
  // itree = 36
  fBoostWeights.push_back(0.0348989);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.507349,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.474517,-99) , 
8, 1.61575, 0, 0, 0.489074,-99)    );
  // itree = 37
  fBoostWeights.push_back(0.0285095);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503521,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.468848,-99) , 
8, 1.7625, 1, 0, 0.490056,-99)    );
  // itree = 38
  fBoostWeights.push_back(0.0369609);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.509577,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476484,-99) , 
4, 2.11884, 1, 0, 0.48846,-99)    );
  // itree = 39
  fBoostWeights.push_back(0.0271771);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503834,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.471862,-99) , 
11, 233.887, 0, 0, 0.491006,-99)    );
  // itree = 40
  fBoostWeights.push_back(0.0265813);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502692,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.473901,-99) , 
2, 1.64463, 0, 0, 0.48966,-99)    );
  // itree = 41
  fBoostWeights.push_back(0.0262444);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502023,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.474348,-99) , 
7, 1.34588, 0, 0, 0.489027,-99)    );
  // itree = 42
  fBoostWeights.push_back(0.0227564);
  fForest.push_back( 
NN(
0, 
0, 
-1, 171.235, 1, -1, 0.488624,-99)    );
  // itree = 43
  fBoostWeights.push_back(0.024546);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.505342,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.47661,-99) , 
1, 171.235, 1, 0, 0.494311,-99)    );
  // itree = 44
  fBoostWeights.push_back(0.029643);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.506381,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.472345,-99) , 
1, 99.6058, 0, 0, 0.492884,-99)    );
  // itree = 45
  fBoostWeights.push_back(0.0284241);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.506265,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.479465,-99) , 
8, 1.61575, 0, 0, 0.491347,-99)    );
  // itree = 46
  fBoostWeights.push_back(0.0309812);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.506236,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.469918,-99) , 
8, 1.7625, 1, 0, 0.492149,-99)    );
  // itree = 47
  fBoostWeights.push_back(0.0261863);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503442,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476943,-99) , 
0, 269.655, 1, 0, 0.490407,-99)    );
  // itree = 48
  fBoostWeights.push_back(0.0235028);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501694,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.472903,-99) , 
1, 99.6058, 0, 0, 0.490298,-99)    );
  // itree = 49
  fBoostWeights.push_back(0.0285537);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503775,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.477382,-99) , 
8, 1.61575, 0, 0, 0.489071,-99)    );
  // itree = 50
  fBoostWeights.push_back(0.0279912);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503164,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.468921,-99) , 
8, 1.7625, 1, 0, 0.489882,-99)    );
  // itree = 51
  fBoostWeights.push_back(0.0277497);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502431,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476832,-99) , 
3, 2.69084, 0, 0, 0.488308,-99)    );
  // itree = 52
  fBoostWeights.push_back(0.0364854);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.508561,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.474652,-99) , 
3, 2.39256, 1, 0, 0.489021,-99)    );
  // itree = 53
  fBoostWeights.push_back(0.0236744);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501739,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.469882,-99) , 
3, 2.24342, 0, 0, 0.490406,-99)    );
  // itree = 54
  fBoostWeights.push_back(0.0290363);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503629,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.47686,-99) , 
8, 1.61575, 0, 0, 0.488695,-99)    );
  // itree = 55
  fBoostWeights.push_back(0.0253936);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501817,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.470119,-99) , 
8, 1.7625, 1, 0, 0.489531,-99)    );
  // itree = 56
  fBoostWeights.push_back(0.0279105);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.50242,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.477563,-99) , 
3, 2.39256, 1, 0, 0.4881,-99)    );
  // itree = 57
  fBoostWeights.push_back(0.0244845);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501086,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.467617,-99) , 
3, 2.24342, 0, 0, 0.489158,-99)    );
  // itree = 58
  fBoostWeights.push_back(0.0252136);
  fForest.push_back( 
NN(
0, 
0, 
-1, 100.407, 1, -1, 0.487396,-99)    );
  // itree = 59
  fBoostWeights.push_back(0.0227922);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.504748,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.480924,-99) , 
6, 100.407, 1, 0, 0.493697,-99)    );
  // itree = 60
  fBoostWeights.push_back(0.0248626);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.506466,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.482854,-99) , 
8, 1.61575, 0, 0, 0.493283,-99)    );
  // itree = 61
  fBoostWeights.push_back(0.0259011);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.505677,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.475558,-99) , 
8, 1.7625, 1, 0, 0.494007,-99)    );
  // itree = 62
  fBoostWeights.push_back(0.0203085);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502236,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.477788,-99) , 
9, 1.79171, 1, 0, 0.492547,-99)    );
  // itree = 63
  fBoostWeights.push_back(0.0206027);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501455,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.475513,-99) , 
9, 1.19479, 0, 0, 0.491492,-99)    );
  // itree = 64
  fBoostWeights.push_back(0.0245697);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502921,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.480322,-99) , 
8, 1.61575, 0, 0, 0.490296,-99)    );
  // itree = 65
  fBoostWeights.push_back(0.0229166);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502016,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.473625,-99) , 
8, 1.7625, 1, 0, 0.491013,-99)    );
  // itree = 66
  fBoostWeights.push_back(0.0248872);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502407,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.479372,-99) , 
3, 2.69084, 0, 0, 0.489722,-99)    );
  // itree = 67
  fBoostWeights.push_back(0.0327284);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.507911,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.477424,-99) , 
3, 2.39256, 1, 0, 0.49035,-99)    );
  // itree = 68
  fBoostWeights.push_back(0.0208467);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501562,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.473572,-99) , 
3, 2.24342, 0, 0, 0.491589,-99)    );
  // itree = 69
  fBoostWeights.push_back(0.0253834);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503278,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.480389,-99) , 
3, 2.39256, 1, 0, 0.490089,-99)    );
   return;
};
 
// Clean up
inline void ReadBDT::Clear() 
{
   for (unsigned int itree=0; itree<fForest.size(); itree++) { 
      delete fForest[itree]; 
   }
}
   inline double ReadBDT::GetMvaValue( const std::vector<double>& inputValues ) const
   {
      // classifier response value
      double retval = 0;

      // classifier response, sanity check first
      if (!IsStatusClean()) {
         std::cout << "Problem in class \"" << fClassName << "\": cannot return classifier response"
                   << " because status is dirty" << std::endl;
         retval = 0;
      }
      else {
         if (IsNormalised()) {
            // normalise variables
            std::vector<double> iV;
            int ivar = 0;
            for (std::vector<double>::const_iterator varIt = inputValues.begin();
                 varIt != inputValues.end(); varIt++, ivar++) {
               iV.push_back(NormVariable( *varIt, fVmin[ivar], fVmax[ivar] ));
            }
            retval = GetMvaValue__( iV );
         }
         else {
            retval = GetMvaValue__( inputValues );
         }
      }

      return retval;
   }
