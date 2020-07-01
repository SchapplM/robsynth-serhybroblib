% Calculate Gravitation load on the joints for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisDE1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:15
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbarprisDE1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbarprisDE1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:14:56
% EndTime: 2020-06-27 17:14:58
% DurationCPUTime: 0.22s
% Computational Cost: add. (262->43), mult. (217->55), div. (14->4), fcn. (28->2), ass. (0->23)
t64 = (qJ(1) ^ 2);
t66 = (pkin(3) ^ 2);
t85 = 2 * pkin(3) * qJ(1) + t64 + t66;
t59 = qJ(1) + pkin(3);
t77 = -pkin(2) + t59;
t78 = -pkin(2) - t59;
t71 = sqrt(-((pkin(1) + t78) * (pkin(1) + t77) * (pkin(1) - t77) * (pkin(1) - t78)));
t84 = g(1) * t71;
t83 = g(2) * t71;
t67 = pkin(2) ^ 2;
t69 = pkin(1) ^ 2;
t80 = -t67 + t69;
t79 = t67 + t69;
t76 = t69 - t85;
t73 = t66 ^ 2;
t65 = pkin(3) * t66;
t63 = qJ(1) * t64;
t62 = pkin(1) - pkin(2);
t61 = pkin(1) + pkin(2);
t60 = -3 * t64;
t51 = t80 + t85;
t50 = t67 + t76;
t1 = [(((g(2) * t50 + t84) * MDP(8) - (g(1) * t50 - t83) * MDP(9)) / pkin(2) * t59 + (-(-(t65 + 4 * t66 * qJ(1) + (5 * t64 + t80) * pkin(3) + 2 * t63) * t84 + g(2) * (6 * t73 * qJ(1) + 2 * (7 * t64 - t79) * t65 + (-6 * qJ(1) * t79 + 16 * t63) * t66 - 2 * t63 * (-t64 + t79) + (t73 + (t62 ^ 2 + t60) * (t61 ^ 2 + t60)) * pkin(3))) * MDP(6) / 0.2e1 + ((MDP(2) / 0.2e1 + MDP(5) / 0.2e1) * (g(2) * t51 - t84) + (-MDP(3) / 0.2e1 + MDP(4) / 0.2e1) * (g(1) * t51 + t83)) * (-t67 + t76)) / (t59 ^ 2)) / pkin(1) / t71;];
taug = t1;
