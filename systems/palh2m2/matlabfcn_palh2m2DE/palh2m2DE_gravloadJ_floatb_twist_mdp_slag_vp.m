% Calculate Gravitation load on the joints for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m2DE_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:56
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m2DE_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:56:46
% EndTime: 2020-06-30 17:56:47
% DurationCPUTime: 0.04s
% Computational Cost: add. (40->19), mult. (80->28), div. (0->0), fcn. (64->8), ass. (0->12)
t22 = sin(qJ(1));
t26 = cos(qJ(1));
t16 = g(1) * t22 - g(2) * t26;
t17 = g(1) * t26 + g(2) * t22;
t19 = sin(qJ(4));
t23 = cos(qJ(4));
t27 = (t16 * t23 + t19 * t17) * MDP(21) + (-t19 * t16 + t23 * t17) * MDP(22);
t25 = cos(qJ(2));
t24 = cos(qJ(3));
t21 = sin(qJ(2));
t20 = sin(qJ(3));
t1 = [t17 * MDP(3) + (-t21 * MDP(10) + t24 * MDP(17) - t20 * MDP(18) + t25 * MDP(9) + MDP(11) + MDP(19) + MDP(2)) * t16 + t27; (-t25 * g(3) + t17 * t21) * MDP(9) + (g(3) * t21 + t17 * t25) * MDP(10); (-t24 * g(3) + t17 * t20) * MDP(17) + (g(3) * t20 + t17 * t24) * MDP(18); t27;];
taug = t1;
