% Calculate joint inertia matrix for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m2DE_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2DE_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'palh2m2DE_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:47
% EndTime: 2020-06-30 17:56:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (50->26), mult. (89->41), div. (0->0), fcn. (55->6), ass. (0->16)
t31 = sin(qJ(4));
t34 = cos(qJ(4));
t48 = MDP(21) * t31 + MDP(22) * t34;
t35 = cos(qJ(3));
t47 = 0.2e1 * t35;
t36 = cos(qJ(2));
t46 = 0.2e1 * t36;
t32 = sin(qJ(3));
t43 = t48 * pkin(5) * t32;
t33 = sin(qJ(2));
t42 = t48 * pkin(4) * t33;
t39 = t36 * pkin(4) + pkin(1);
t28 = pkin(2) + t39;
t38 = t35 * pkin(5) + t28;
t37 = (t34 * MDP(21) - t31 * MDP(22)) * (pkin(3) + t38);
t1 = [pkin(1) * MDP(9) * t46 + t28 * MDP(17) * t47 + MDP(1) + MDP(20) + 0.2e1 * t37 + 0.2e1 * t39 * MDP(11) + 0.2e1 * t38 * MDP(19) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t33 + MDP(5) * t46) * t33 + (MDP(12) * t32 + MDP(13) * t47 - 0.2e1 * t28 * MDP(18)) * t32; t33 * MDP(6) + t36 * MDP(7) + t42; MDP(8); t32 * MDP(14) + t35 * MDP(15) + t43; ((t33 * t32 + t36 * t35) * MDP(17) + (-t36 * t32 + t35 * t33) * MDP(18)) * pkin(4); MDP(16); MDP(20) + t37; t42; t43; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
