% Calculate Gravitation load on the joints for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m2OL_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m2OL_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'palh2m2OL_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:06:39
% EndTime: 2020-06-30 18:06:39
% DurationCPUTime: 0.12s
% Computational Cost: add. (362->44), mult. (278->62), div. (0->0), fcn. (244->12), ass. (0->29)
t56 = qJ(2) + qJ(3);
t55 = qJ(4) + t56;
t52 = qJ(5) + t55;
t48 = sin(t52);
t74 = g(3) * t48;
t57 = sin(qJ(6));
t59 = sin(qJ(1));
t72 = t59 * t57;
t60 = cos(qJ(6));
t71 = t59 * t60;
t62 = cos(qJ(1));
t70 = t62 * t57;
t69 = t62 * t60;
t49 = cos(t52);
t67 = g(1) * t62 + g(2) * t59;
t68 = (t67 * t49 + t74) * MDP(31) + (t60 * MDP(37) - t57 * MDP(38) + MDP(30)) * (-g(3) * t49 + t67 * t48);
t50 = sin(t55);
t51 = cos(t55);
t65 = (-g(3) * t51 + t67 * t50) * MDP(23) + (g(3) * t50 + t67 * t51) * MDP(24) + t68;
t53 = sin(t56);
t54 = cos(t56);
t64 = (-g(3) * t54 + t67 * t53) * MDP(16) + (g(3) * t53 + t67 * t54) * MDP(17) + t65;
t61 = cos(qJ(2));
t58 = sin(qJ(2));
t47 = t49 * t69 - t72;
t46 = -t49 * t70 - t71;
t45 = -t49 * t71 - t70;
t44 = t49 * t72 - t69;
t1 = [t67 * MDP(3) + (-g(1) * t45 - g(2) * t47) * MDP(37) + (-g(1) * t44 - g(2) * t46) * MDP(38) + (-t58 * MDP(10) + MDP(16) * t54 - MDP(17) * t53 + MDP(23) * t51 - MDP(24) * t50 + MDP(30) * t49 - MDP(31) * t48 + t61 * MDP(9) + MDP(2)) * (g(1) * t59 - g(2) * t62); (-g(3) * t61 + t67 * t58) * MDP(9) + (g(3) * t58 + t67 * t61) * MDP(10) + t64; t64; t65; t68; (-g(1) * t46 + g(2) * t44 + t57 * t74) * MDP(37) + (g(1) * t47 - g(2) * t45 + t60 * t74) * MDP(38);];
taug = t1;
