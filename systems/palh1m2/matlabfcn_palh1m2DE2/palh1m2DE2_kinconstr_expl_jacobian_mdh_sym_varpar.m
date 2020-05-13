% Jacobian of explicit kinematic constraints of
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% W [16x4]
%  Derivative of the joint coordinates w.r.t minimal coordinates
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function W = palh1m2DE2_kinconstr_expl_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_kinconstr_expl_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_kinconstr_expl_jacobian_mdh_sym_varpar: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:24
% EndTime: 2020-05-02 20:57:25
% DurationCPUTime: 0.23s
% Computational Cost: add. (750->0), mult. (1284->0), div. (70->0), fcn. (1920->0), ass. (0->1)
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, -1, -1, 0; 0, 0, 0, 1; 0, 1, 0, 0; 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, -1, 0; 0, 1, 1, 0; 0, 1, 0, 0; 0, 0, -1, 0; 0, 1, 1, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
W = t1;
