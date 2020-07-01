% Calculate joint inertia matrix for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m1DE_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1DE_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'palh2m1DE_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:36
% EndTime: 2020-06-30 17:39:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (72->35), mult. (140->52), div. (0->0), fcn. (113->6), ass. (0->21)
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t66 = MDP(20) * t46 + MDP(21) * t49;
t50 = cos(qJ(3));
t44 = t50 * pkin(3) + pkin(2);
t51 = cos(qJ(2));
t47 = sin(qJ(3));
t48 = sin(qJ(2));
t63 = t48 * t47;
t53 = -pkin(3) * t63 + t44 * t51 + pkin(1);
t65 = (t49 * MDP(20) - t46 * MDP(21)) * (pkin(4) + t53);
t62 = t48 * t51;
t61 = t51 * t47;
t60 = t66 * (pkin(3) * t61 + t48 * t44);
t38 = -t50 * t48 - t61;
t59 = t66 * pkin(3) * t38;
t40 = -t51 * t50 + t63;
t58 = t38 * MDP(13) + t40 * MDP(14);
t52 = (MDP(16) * t50 - MDP(17) * t47) * pkin(2);
t45 = t51 * pkin(2) + pkin(1);
t1 = [MDP(1) + t48 ^ 2 * MDP(4) + 0.2e1 * MDP(5) * t62 + (0.4e1 * (t50 ^ 2 - 0.1e1 / 0.2e1) * t62 + (0.4e1 * t51 ^ 2 - 0.2e1) * t50 * t47) * MDP(12) - 0.2e1 * t45 * t40 * MDP(16) + 0.2e1 * t53 * MDP(18) + MDP(19) + (MDP(11) * t38 + 0.2e1 * t45 * MDP(17)) * t38 + 0.2e1 * t65 + 0.2e1 * (-t48 * MDP(10) + t51 * MDP(9)) * pkin(1); -t48 * MDP(6) - t51 * MDP(7) + t58 + t60; MDP(15) + MDP(8) + 0.2e1 * t52; t58 - t59; MDP(15) + t52; MDP(15); MDP(19) + t65; t60; -t59; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
