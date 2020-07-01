% Calculate Gravitation load on the joints for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m1OL_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:46
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m1OL_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:45:12
% EndTime: 2020-06-30 17:45:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (199->37), mult. (202->54), div. (0->0), fcn. (180->10), ass. (0->25)
t46 = qJ(2) + qJ(3);
t45 = qJ(4) + t46;
t41 = sin(t45);
t63 = g(3) * t41;
t47 = sin(qJ(5));
t49 = sin(qJ(1));
t61 = t49 * t47;
t50 = cos(qJ(5));
t60 = t49 * t50;
t52 = cos(qJ(1));
t59 = t52 * t47;
t58 = t52 * t50;
t42 = cos(t45);
t56 = g(1) * t52 + g(2) * t49;
t57 = (t42 * t56 - t63) * MDP(24) + (t50 * MDP(30) - t47 * MDP(31) + MDP(23)) * (g(3) * t42 + t41 * t56);
t43 = sin(t46);
t44 = cos(t46);
t54 = (g(3) * t44 + t43 * t56) * MDP(16) + (-g(3) * t43 + t44 * t56) * MDP(17) + t57;
t51 = cos(qJ(2));
t48 = sin(qJ(2));
t40 = t42 * t58 - t61;
t39 = -t42 * t59 - t60;
t38 = -t42 * t60 - t59;
t37 = t42 * t61 - t58;
t1 = [t56 * MDP(3) + (-g(1) * t38 - g(2) * t40) * MDP(30) + (-g(1) * t37 - g(2) * t39) * MDP(31) + (-MDP(10) * t48 + MDP(16) * t44 - MDP(17) * t43 + MDP(23) * t42 - MDP(24) * t41 + MDP(9) * t51 + MDP(2)) * (g(1) * t49 - g(2) * t52); (g(3) * t51 + t48 * t56) * MDP(9) + (-g(3) * t48 + t51 * t56) * MDP(10) + t54; t54; t57; (-g(1) * t39 + g(2) * t37 - t47 * t63) * MDP(30) + (g(1) * t40 - g(2) * t38 - t50 * t63) * MDP(31);];
taug = t1;
