% Calculate joint inertia matrix for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see fourbarprisDE1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [1x1]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:15
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbarprisDE1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_inertiaJ_mdp_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_inertiaJ_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'fourbarprisDE1_inertiaJ_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:14:55
% EndTime: 2020-06-27 17:14:55
% DurationCPUTime: 0.04s
% Computational Cost: add. (87->17), mult. (47->20), div. (20->3), fcn. (2->2), ass. (0->8)
t36 = (qJ(1) + pkin(3));
t39 = -pkin(2) + t36;
t40 = -pkin(2) - t36;
t45 = ((pkin(1) - t39) * (pkin(1) + t40) * (pkin(1) - t40) * (pkin(1) + t39));
t38 = (t36 ^ 2);
t37 = (qJ(1) ^ 2);
t25 = pkin(1) ^ 2 - pkin(2) ^ 2 - t37 - (2 * qJ(1) + pkin(3)) * pkin(3);
t1 = [2 * t25 / t36 * (-t45) ^ (-0.1e1 / 0.2e1) * MDP(4) + MDP(6) + (-4 * MDP(7) * t38 + (-2 * MDP(5) * qJ(1) - t37 * MDP(6) - MDP(1)) / t38 * t25 ^ 2) / t45;];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t1(1);];
Mq = res;
