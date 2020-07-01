% Calculate Gravitation load on the joints for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1TE_convert_par2_MPV_fixb.m
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:21
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1TE_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbar1TE_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:21:05
% EndTime: 2020-06-26 17:21:06
% DurationCPUTime: 0.22s
% Computational Cost: add. (323->43), mult. (498->58), div. (16->5), fcn. (140->4), ass. (0->31)
t70 = cos(qJ(1));
t88 = pkin(2) * t70;
t82 = (-0.2e1 * t88 + pkin(1)) * pkin(1);
t66 = pkin(2) ^ 2 + t82;
t64 = 0.1e1 / t66;
t92 = t64 / 0.2e1;
t91 = -pkin(3) - pkin(4);
t90 = -pkin(3) + pkin(4);
t69 = sin(qJ(1));
t89 = pkin(2) * t69;
t56 = (pkin(2) - t91) * (pkin(2) + t91) + t82;
t57 = (pkin(2) - t90) * (pkin(2) + t90) + t82;
t78 = sqrt(-t56 * t57);
t80 = pkin(1) * t89;
t87 = (-t56 - t57) * t80 / t78;
t62 = -g(1) * t89 + g(2) * t88;
t58 = g(2) * pkin(1) - t62;
t86 = t58 * t64;
t79 = g(1) * t70 + g(2) * t69;
t63 = t79 * pkin(2);
t59 = -g(1) * pkin(1) + t63;
t85 = t59 * t64;
t84 = -t58 * t87 - t63 * t78;
t83 = -t59 * t87 - t62 * t78;
t81 = pkin(3) ^ 2 - pkin(4) ^ 2;
t65 = 0.1e1 / t66 ^ 2;
t61 = t66 - t81;
t60 = t66 + t81;
t50 = t59 * t78;
t49 = t58 * t78;
t1 = [(g(1) * t69 - g(2) * t70) * MDP(2) + t79 * MDP(3) + (((t60 * t62 + t84) * t92 + (t85 - (t60 * t59 - t49) * t65) * t80) * MDP(5) + ((-t60 * t63 + t83) * t92 + (-t86 - (-t60 * t58 - t50) * t65) * t80) * MDP(6)) / pkin(3) + (((-t61 * t62 + t84) * t92 + (-t85 - (-t61 * t59 - t49) * t65) * t80) * MDP(8) + ((t61 * t63 + t83) * t92 + (t86 - (t61 * t58 - t50) * t65) * t80) * MDP(9)) / pkin(4);];
taug = t1;
