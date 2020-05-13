% Explicit kinematic constraints of
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% jv [5x1]
%   Joint variables (rotation around z or translation in z-direction according to MDH)
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function jv = fourbarprisDE2_kinconstr_expl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_kinconstr_expl_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_kinconstr_expl_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:43:51
% EndTime: 2020-05-07 09:43:51
% DurationCPUTime: 0.08s
% Computational Cost: add. (99->16), mult. (66->24), div. (20->4), fcn. (11->5), ass. (0->17)
t17 = -qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t22 = (pkin(1) ^ 2 - pkin(2) ^ 2);
t3 = t17 + t22;
t25 = -t3 / 0.4e1;
t14 = 1 / pkin(2);
t10 = qJ(1) + pkin(3);
t8 = 1 / t10;
t24 = t14 * t8;
t16 = 1 / pkin(1);
t23 = t16 * t8;
t21 = t14 * t16 / t10 ^ 2;
t20 = -pkin(2) - t10;
t19 = -pkin(2) + t10;
t18 = (pkin(1) + t20) * (pkin(1) + t19) * (pkin(1) - t19) * (pkin(1) - t20);
t2 = -t17 + t22;
t1 = sqrt(-t18);
t4 = [atan2(t1 * t23, (t2 * t23)); qJ(1); atan2((t25 + t2 / 0.4e1) * t1 * t21, (t2 * t25 + t18 / 0.4e1) * t21); atan2(t1 * t24, (t3 * t24)); -pi;];
jv = t4(:);
