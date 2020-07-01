% Calculate Gravitation load on the joints for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see palh2m1DE_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m1DE_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:36
% EndTime: 2020-06-30 17:39:36
% DurationCPUTime: 0.04s
% Computational Cost: add. (58->23), mult. (113->36), div. (0->0), fcn. (98->8), ass. (0->15)
t36 = cos(qJ(1));
t30 = sin(qJ(1));
t25 = g(1) * t36 + g(2) * t30;
t28 = sin(qJ(3));
t32 = cos(qJ(3));
t22 = t32 * g(3) + t28 * t25;
t23 = -t28 * g(3) + t25 * t32;
t29 = sin(qJ(2));
t33 = cos(qJ(2));
t35 = (t22 * t33 + t23 * t29) * MDP(16) + (-t29 * t22 + t23 * t33) * MDP(17);
t24 = g(1) * t30 - g(2) * t36;
t27 = sin(qJ(4));
t31 = cos(qJ(4));
t34 = (t24 * t31 + t27 * t25) * MDP(20) + (-t27 * t24 + t31 * t25) * MDP(21);
t1 = [t25 * MDP(3) + (MDP(2) + t33 * MDP(9) - t29 * MDP(10) + (-t29 * t28 + t33 * t32) * MDP(16) - (t33 * t28 + t32 * t29) * MDP(17) + MDP(18)) * t24 + t34; (g(3) * t33 + t25 * t29) * MDP(9) + (-g(3) * t29 + t25 * t33) * MDP(10) + t35; t35; t34;];
taug = t1;
