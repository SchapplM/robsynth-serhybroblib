% Calculate joint inertia matrix for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbar1turnOL_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1turnOL_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'fourbar1turnOL_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:30
% EndTime: 2020-06-27 16:56:30
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->20), mult. (65->35), div. (0->0), fcn. (54->6), ass. (0->14)
t19 = sin(qJ(3));
t20 = sin(qJ(2));
t22 = cos(qJ(3));
t23 = cos(qJ(2));
t16 = t19 * t20 - t22 * t23;
t28 = 0.2e1 * t16;
t21 = cos(qJ(4));
t27 = 0.2e1 * t21;
t26 = pkin(2) * t23;
t17 = -t19 * t23 - t22 * t20;
t25 = t17 * MDP(13) + t16 * MDP(14);
t24 = (-MDP(16) * t22 + MDP(17) * t19) * pkin(2);
t18 = sin(qJ(4));
t1 = [MDP(16) * t26 * t28 + pkin(1) * MDP(23) * t27 + MDP(1) + (MDP(4) * t20 + 0.2e1 * MDP(5) * t23) * t20 + (MDP(18) * t18 + MDP(19) * t27 - 0.2e1 * MDP(24) * pkin(1)) * t18 + (MDP(11) * t17 + MDP(12) * t28 - 0.2e1 * MDP(17) * t26) * t17; t20 * MDP(6) + t23 * MDP(7) + t25; MDP(15) + MDP(8) + 0.2e1 * t24; t25; MDP(15) + t24; MDP(15); t18 * MDP(20) + t21 * MDP(21); 0; 0; MDP(22); 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
